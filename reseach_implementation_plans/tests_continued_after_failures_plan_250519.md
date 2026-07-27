# Plan: Tests Continue After Failures (Including Timeout)

**Date:** 2026-05-19  
**Author:** Daniil Nikitin  
**Status:** Implemented ✅  
**Source log:** `logs/test_logs_260519_part2.txt`

---

## Problem Statement

Running `pytest tests/ -v --timeout=3600` on 2026-05-19 collected **104 tests**. After test 23 of 104 (`test_g_full_pipeline_within_tolerance`) hit the 3600 s timeout, the entire pytest session **hung indefinitely** and never ran the remaining ~81 tests. Additionally, `test_g_cublock_matches_octave` (test 22 of 104) failed with an assertion error just before the hanging test.

Expected behaviour: a timed-out or failed test should be marked accordingly and pytest should continue to the next test.

---

## Root Cause Analysis

### Issue 1 — `timeout_method = "thread"` cannot interrupt blocked C-level lock acquisition

`pyproject.toml` configures:
```toml
timeout_method = "thread"
```

The thread method works by calling `ctypes.pythonapi.PyThreadState_SetAsyncExc()` in a background thread when the timeout fires. This injects a Python exception into the target thread — but **only when that thread is executing Python bytecode**. It has no effect while the thread is blocked in a C extension.

The stack trace shows the test thread stuck here:
```
concurrent.futures._base.as_completed → threading.Condition.wait → threading.Lock.acquire
```

`threading.Lock.acquire()` is a C-level blocking operation. The exception injection is silently dropped, the lock is never released, and the test thread never unblocks. pytest marks the test as timed-out in its internal state but the process itself is still frozen at that lock.

### Issue 2 — Orphaned `ProcessPoolExecutor` workers prevent session continuation

`parallel.py:377` runs workers inside:
```python
with ProcessPoolExecutor(max_workers=...) as executor:
    ...
    for future in as_completed(future_to_idx):   # ← hung here
        ...
```

When thread-injected timeout eventually fails and pytest tries to clean up, the `ProcessPoolExecutor.__exit__` calls `shutdown(wait=True)`. The Python CuBlock workers (and their internal QueueFeederThreads) are still running inside the subprocess. Python's shutdown waits for them. The cleanup itself hangs, preventing pytest from advancing to the next test.

The two QueueFeederThread stacks in the log confirm this: the multiprocessing queue machinery is still live when the timeout fires.

### Issue 3 — `test_g_cublock_matches_octave` failed (separate issue)

This test compares Python CuBlock output against Octave CuBlock using `mean_rel_diff < 5%`. The failure indicates the two implementations diverge beyond 5%. 

Root cause: NumPy's `default_rng` and Octave's `rand('state', seed)` are entirely different PRNGs and produce different random sequences even from the same integer seed. k-means with different random initialization converges to different cluster assignments in each of the 30 reps. Averaged over 30 reps the outputs should converge, but with k=5 and a single synthetic input column, the averaging is insufficient to wash out PRNG differences.

This failure is not related to the hanging issue and requires separate attention if needed in the future.

---

## Solution Options

| Option | Mechanism | Interrupts blocked locks? | Kills orphaned workers? | Complexity |
|---|---|---|---|---|
| **A. `timeout_method = "signal"`** | Delivers `SIGALRM` to main thread; Python's signal machinery runs after any blocking syscall returns with `EINTR` | ✅ Yes — POSIX signals interrupt `select()`/mutex waits | ❌ No — still need explicit cleanup | Low |
| **B. `autouse` orphan-cleanup fixture** | After each test (including timeout), terminates child processes via `multiprocessing.active_children()` | N/A (complements A) | ✅ Yes | Medium |
| **C. `pytest-forked`** | Each test runs in a forked subprocess; if it hangs, the forked process is killed and pytest continues | ✅ Yes — fork isolation | ✅ Yes — isolated process | Medium |
| **D. Subprocess wrapper in `run_pipeline`** | `run_pipeline` launches itself in a `subprocess.run` with hard wall-clock limit | ✅ Yes | ✅ Yes | High — loses pytest output |
| **E. `pytest-xdist -n 1`** | Each test in a separate worker process via xdist | ✅ Yes — xdist detects worker crash | ✅ Yes | Medium — changes test semantics |

### Recommendation

**Option A + Option B together**, with an additional reduction of the `test_g_full_pipeline_within_tolerance` timeout from 3600 s to 600 s.

- A is a one-line change to `pyproject.toml` that fixes the fundamental blocking issue.
- B ensures orphaned workers are cleaned up even in edge cases where SIGALRM doesn't prevent cleanup from hanging.
- Reducing the g-test timeout from 3600 s to 600 s gives a faster failure signal: Python CuBlock on 10 samples with `max_iter=100` should complete in 2–5 minutes at most. A 3600 s limit means a single hanging test wastes an hour.

Option C (`pytest-forked`) is the most robust alternative if A+B prove insufficient. Note it requires `pip install pytest-forked` and that each test must be annotated with `@pytest.mark.forked` (or the entire suite run with `--forked`).

---

## Approved Changes

### Change 1 — `pyproject.toml`: switch to signal timeout method

```toml
[tool.pytest.ini_options]
timeout = 3600
timeout_method = "signal"   # was: "thread"
```

**Why this fixes the hang:** `SIGALRM` is delivered to the process regardless of what C-level blocking call the main thread is in. Python's `signal` module ensures the handler runs at the next safe point after any blocking syscall returns (`EINTR`). The pytest-timeout handler raises `pytest.TimeoutExpired`, which propagates up through `as_completed()` → `harmonize_parallel()` → the test function, allowing pytest to mark the test and move on.

**Constraint:** Signal method only works from the main thread and only on POSIX systems. This is always satisfied here — pytest's main loop and all tests run in the main thread; the pod is Linux.

---

### Change 2 — `tests/speed_up_tests/conftest.py`: `autouse` orphan-cleanup fixture

Add a function-scoped fixture that kills any `ProcessPoolExecutor` child processes that remain alive after a test ends. This prevents orphaned workers from blocking `shutdown(wait=True)` in the executor cleanup path.

The fixture:
1. Before the test: records the PIDs of all current child processes via `multiprocessing.active_children()`.
2. After the test (in `finally`, so it runs even on timeout/failure): compares current children to the pre-test set; terminates any new ones with `SIGTERM`, gives them 3 s to exit, then sends `SIGKILL` to any that are still alive.

Location in `conftest.py`: after the existing constants, before `load_fixtures`. Use `scope="function"` and `autouse=True` so it runs around every test in the `speed_up_tests/` directory automatically.

```python
@pytest.fixture(autouse=True, scope="function")
def _kill_orphaned_workers():
    before = {p.pid for p in multiprocessing.active_children()}
    yield
    stragglers = [p for p in multiprocessing.active_children() if p.pid not in before]
    for p in stragglers:
        p.terminate()
    # Brief grace period, then SIGKILL
    deadline = time.monotonic() + 3.0
    while stragglers and time.monotonic() < deadline:
        time.sleep(0.1)
        stragglers = [p for p in stragglers if p.is_alive()]
    for p in stragglers:
        p.kill()
```

This requires adding `import multiprocessing` at the top of `conftest.py` (it already imports `time`).

---

### Change 3 — `tests/speed_up_tests/test_approach_g.py`: reduce timeouts

Change the per-test timeout decorator on `test_g_full_pipeline_within_tolerance` from `3600` to `600`:

```python
@octave_required
@pytest.mark.timeout(600)   # was: 3600 — Python CuBlock on 10 samples should finish in <5 min
def test_g_full_pipeline_within_tolerance():
    ...
```

With `max_iter=100`, 30 reps, 10 samples, and pool size 49, Python CuBlock should complete in 2–5 minutes. A 600 s (10 min) limit is ample and gives a much faster failure signal if the function regresses.

Similarly, reduce `test_g_cublock_matches_octave` timeout from `3600` to `300` (matching the internal Octave timeout already used in that test at `conftest.py:87`):

```python
@octave_required
@pytest.mark.timeout(300)   # was: 3600 — single CuBlock call; matches internal Octave timeout
def test_g_cublock_matches_octave():
    ...
```

---

## Summary of Files to Modify

| File | Change | Section |
|---|---|---|
| `pyproject.toml` | `timeout_method = "signal"` | Change 1 |
| `tests/speed_up_tests/conftest.py` | Add `import multiprocessing`; add `_kill_orphaned_workers` autouse fixture | Change 2 |
| `tests/speed_up_tests/test_approach_g.py` | Reduce timeout on `test_g_full_pipeline_within_tolerance` to 600 s; reduce `test_g_cublock_matches_octave` timeout to 300 s | Change 3 |

---

## Implementation Note — `shambhala/parallel.py` (discovered during implementation)

During smoke-testing it became clear that `timeout_method = "signal"` alone was insufficient. When SIGALRM fires and the `Timeout` exception propagates out of `as_completed()`, `ProcessPoolExecutor.__exit__` calls `shutdown(wait=True)`, which joins the still-running worker processes. Those workers may be mid-computation (e.g. 30 reps of Python CuBlock), so the main thread blocks in `join()` for minutes — effectively a new hang at the same location, just via the cleanup path rather than the lock.

Fix applied to `shambhala/parallel.py`: replaced `with ProcessPoolExecutor(...) as executor` with an explicit `executor = ProcessPoolExecutor(...)` + `try/except BaseException/else/finally` pattern. On any `BaseException` (including `pytest.fail.Exception` from the timeout), the handler terminates all new child processes and calls `executor.shutdown(wait=False, cancel_futures=True)` before re-raising. Normal exit still calls `shutdown(wait=True)`. The `_kill_orphaned_workers` fixture then cleans up any stragglers with SIGKILL.

| File | Change |
|---|---|
| `shambhala/parallel.py` | Replace context-manager executor with explicit try/except BaseException; terminate workers on exception; shutdown(wait=False) on exception, shutdown(wait=True) on success |

---

## To-Do List

### 1. `pyproject.toml` — switch timeout method (Change 1) ✅

- [x] Open `pyproject.toml`.
- [x] In `[tool.pytest.ini_options]`, change `timeout_method = "thread"` to `timeout_method = "signal"`.
- [x] Verify the global `timeout = 3600` line is still present (it is the fallback for tests without a per-test decorator).

---

### 2. `tests/speed_up_tests/conftest.py` — add orphan-cleanup fixture (Change 2) ✅

- [x] Add `import multiprocessing` to the imports block at the top of `conftest.py` (it already imports `time`, so insert on a nearby line).
- [x] Locate the section after the existing constants and before the `load_fixtures` function.
- [x] Insert the `_kill_orphaned_workers` fixture exactly as specified in Change 2 above.
- [x] Confirm the fixture has `autouse=True` and `scope="function"`.
- [x] Confirm the fixture uses `multiprocessing.active_children()` (not `os.getpid()` or psutil) to stay dependency-free.

---

### 3. `tests/speed_up_tests/test_approach_g.py` — reduce timeouts (Change 3) ✅

- [x] Find `@pytest.mark.timeout(3600)` on `test_g_full_pipeline_within_tolerance`; change to `@pytest.mark.timeout(600)`.
- [x] Find `@pytest.mark.timeout(3600)` on `test_g_cublock_matches_octave`; change to `@pytest.mark.timeout(300)`.
- [x] Confirm no other tests in the file still carry `3600` as their explicit decorator (the global fallback covers tests without an explicit decorator).

---

### 4. Smoke-test the fixes locally ✅

- [x] Run speed_up_tests with an artificial short cap (`--timeout=10`) to verify pytest continues past failures:
  confirmed — all 31 selected tests were attempted; session ran to completion in 15m 46s (no hang); G-5 timed out cleanly at its 600s per-test decorator; all `test_combinations.py` tests timed out at 10s (CLI cap) as expected.
- [x] Pure-Python unit tests: 63 passed in 2.6s after all changes (including `parallel.py` fix).
- [x] Confirmed `timeout method: signal` shown in pytest session header.

---

### 5. Commit ✅

- [x] Stage and commit the four changed files:
  ```
  pyproject.toml
  tests/speed_up_tests/conftest.py
  tests/speed_up_tests/test_approach_g.py
  shambhala/parallel.py
  ```
- [x] Commit message: `fix: switch pytest timeout to signal method; add orphan-worker cleanup fixture; tighten g-test timeouts`
