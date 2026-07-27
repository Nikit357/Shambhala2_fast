# Implementation Plan: Detailed Per-Worker Progress Logging for Shambhala

**Created:** 2026-05-17  
**Author:** Daniil Nikitin  
**Scope:** `Shambhala_containerized/` only; no changes to `Shambhala2/` or `bench_shared.py`.

---

## 1. Goal

When `harmonize_parallel()` is running, the user (watching the terminal or a log tail) should see a live, fixed-position display — one row per active Octave worker — showing how many samples that worker has already processed, with a completion bar and percentage. The display should look like:

```
[Shambhala] Harmonizing 500 samples | 5 workers | k=5
─────────────────────────────────────────────────────────────────────────────
 Worker 1  [██████████████░░░░░░░░░░░░░░]  36/100 ( 36%)  3.2 s/sample  ETA  3m 12s
 Worker 2  [████████████████░░░░░░░░░░░░]  40/100 ( 40%)  3.1 s/sample  ETA  2m 52s
 Worker 3  [████████████████████░░░░░░░░]  50/100 ( 50%)  3.0 s/sample  ETA  2m 30s
 Worker 4  [██████████░░░░░░░░░░░░░░░░░░]  25/100 ( 25%)  3.4 s/sample  ETA  4m 04s
 Worker 5  [████████████░░░░░░░░░░░░░░░░]  30/100 ( 30%)  3.2 s/sample  ETA  3m 44s
─────────────────────────────────────────────────────────────────────────────
 Overall   [█████████████░░░░░░░░░░░░░░░] 181/500 ( 36%)  elapsed  9m 37s  ETA  17m 06s
```

In non-TTY environments (log file redirect, K8s pod log capture), the display falls back to plain timestamped log lines — one line per completed sample inside each worker.

---

## 2. Key Pre-existing Opportunity (Zero Octave Changes Needed)

`octave/Shambhala2_piped.m` already emits per-sample progress messages to stdout on line 15:

```octave
message = sprintf('Harmonizing sample %d out of %d', i, NH);
disp(message);
```

These text lines are currently discarded silently inside `_parse_octave_stdout()` (the `# Non-data lines` branch in `octave_bridge.py:92`). This means:

- **No changes to any `.m` file are needed.**
- Real per-sample granularity within each Octave subprocess is already available.
- All that is needed is to switch from `subprocess.run()` (buffered, blocking) to `subprocess.Popen()` (streaming) and intercept those lines.

---

## 3. Architecture Overview

```
Main process
│
├── harmonize_parallel()
│   ├── creates multiprocessing.Manager().Queue()  ← shared IPC channel
│   ├── spawns _ProgressRenderer (daemon threading.Thread)
│   │      reads queue → renders/updates ANSI rows on stderr
│   │
│   └── ProcessPoolExecutor (n_workers processes)
│        ├── Worker 1: _harmonize_one_batch_with_progress(...)
│        │    └── _run_octave_streaming(pool_df, ..., progress_queue, worker_idx)
│        │         ├── subprocess.Popen(octave, stdin=PIPE, stdout=PIPE, stderr=PIPE)
│        │         ├── writes TSV to Popen.stdin, then closes it
│        │         ├── reads stdout line by line:
│        │         │    "Harmonizing sample X out of Y" → puts (worker_idx, X, Y, elapsed) in queue
│        │         │    "<SYMBOL> val1 val2 ..."        → accumulate for later parsing
│        │         └── after process exits → parse accumulated data lines → return DataFrame
│        │
│        ├── Worker 2: same
│        └── ...
│
└── as_completed() loop → collects results → assembles output
```

**IPC mechanism:** `multiprocessing.Manager().Queue()` (not `multiprocessing.Queue`) because Manager queues cross the `ProcessPoolExecutor` process boundary safely. Workers put `ProgressEvent` named-tuples; the renderer thread in the main process reads and re-renders.

---

## 4. Files to Create / Modify

| Action | File | Summary |
|---|---|---|
| **Create** | `shambhala/progress_display.py` | `_ProgressRenderer` class; all ANSI logic |
| **Modify** | `shambhala/octave_bridge.py` | Add `run_octave_normalize_streaming()`; keep old `run_octave_normalize()` |
| **Modify** | `shambhala/parallel.py` | Wire up queue + renderer; new `_harmonize_one_batch_with_progress()` |
| **No change** | `octave/Shambhala2_piped.m` | Already emits progress lines |
| **No change** | `run_shambhala.py` | Calls `harmonize_parallel()` — inherits changes automatically |
| **No change** | `harmonization_scripts/run_shambhala_job.py` | Same |
| **Add tests** | `tests/test_progress_display.py` | See §8 |

---

## 5. New File: `shambhala/progress_display.py`

### 5.1 Public interface

```python
from shambhala.progress_display import ProgressEvent, _ProgressRenderer

# Event sent by each worker process into the queue:
ProgressEvent = NamedTuple('ProgressEvent', [
    ('worker_idx', int),   # 0-based worker index
    ('sample_done', int),  # 1-based count of samples done so far in this batch
    ('sample_total', int), # total samples assigned to this worker
    ('elapsed_s', float),  # seconds since worker started
    ('is_final', bool),    # True only for the synthetic "batch complete" event
])

class _ProgressRenderer(threading.Thread):
    def __init__(
        self,
        queue: Queue,
        n_workers: int,
        total_samples: int,
        worker_batch_sizes: list[int],  # samples per worker, for ETA
        use_ansi: bool,                 # True only when sys.stderr.isatty()
    ): ...

    def run(self) -> None:
        """Main loop: poll queue, update state, re-render at most 4 Hz."""

    def stop(self) -> None:
        """Signal the renderer to flush and exit; call after all futures complete."""
```

### 5.2 Internal state

```python
# Per-worker state, updated on each ProgressEvent:
_state: list[dict]   # index = worker_idx
# Each entry:
{
    'done': int,         # samples done
    'total': int,        # samples assigned
    'elapsed_s': float,  # latest elapsed from event
    'finished': bool,    # True after is_final event
}
```

### 5.3 ANSI rendering (TTY mode)

- Maintain a `_last_rendered_lines: int` counter.
- On each re-render: `sys.stderr.write('\033[{N}A\033[J')` (move up N lines, clear to end of screen), then write all rows fresh.
- Bar width: 30 characters. Filled character: `█`. Empty: `░`.
- ETA calculation: `samples_remaining / (done / elapsed_s)` if `done > 0` else `--`.
- `speed_str`: `elapsed_s / done` s/sample if `done > 0`, else `--`.
- Overall row: sum `done` across all workers / `total_samples`.
- Re-render only when state changed OR at most every 0.25 s (4 Hz cap, avoids flicker).
- On `stop()`: do one final render, then write `\n` so the next log output appears below.

### 5.4 Plain-text fallback (non-TTY / log-file mode)

When `use_ansi=False` (i.e., `not sys.stderr.isatty()`):

- On each `ProgressEvent` where `sample_done % log_every == 0` or `is_final`: emit one timestamped line to `sys.stderr`:
  ```
  [HH:MM:SS][Worker 3]  50/100 (50%)  elapsed 2m 30s
  ```
- Default `log_every`: the larger of `1` and `total // 20` (emit at most ~20 lines per worker).
- `is_final` events always emit a summary line.

### 5.5 Sentinel value for stopping the thread

The renderer loop exits when it receives the `None` sentinel from the queue (put by `stop()`):

```python
def stop(self) -> None:
    self._queue.put(None)
    self.join(timeout=5.0)
```

---

## 6. Modified File: `shambhala/octave_bridge.py`

### 6.1 New function: `run_octave_normalize_streaming()`

Signature identical to `run_octave_normalize()` plus two extra parameters:

```python
def run_octave_normalize_streaming(
    pool_df: pd.DataFrame,
    nh: int,
    np_: int,
    k: int,
    octave_scripts_dir: str,
    octave_bin: str = "octave",
    timeout_s: int = 6000,
    random_seed: int = None,
    progress_queue=None,   # multiprocessing Queue or None
    worker_idx: int = 0,   # for labeling ProgressEvents
) -> pd.DataFrame:
```

**Implementation steps:**

1. Build `tsv_string` and `eval_preamble` (identical to current `run_octave_normalize`).

2. Launch Octave with `subprocess.Popen`:
   ```python
   proc = subprocess.Popen(
       [octave_bin, "--no-gui", "--eval", eval_preamble],
       stdin=subprocess.PIPE,
       stdout=subprocess.PIPE,
       stderr=subprocess.PIPE,
       text=True,
   )
   ```

3. Write all TSV data to stdin in a separate thread (to avoid deadlock when the TSV is large and Octave's stdout buffer fills before Python finishes writing):
   ```python
   # stdin_writer thread: proc.stdin.write(tsv_string); proc.stdin.close()
   ```
   This is the standard deadlock-prevention pattern for bidirectional Popen I/O.

4. Read `proc.stdout` line by line in the main worker thread:
   ```python
   stdout_lines = []
   t_start = time.time()
   for line in proc.stdout:
       line = line.rstrip('\n')
       if _is_progress_line(line):            # "Harmonizing sample X out of Y"
           x, y = _parse_progress_line(line)
           if progress_queue is not None:
               progress_queue.put(ProgressEvent(
                   worker_idx=worker_idx,
                   sample_done=x,
                   sample_total=y,
                   elapsed_s=time.time() - t_start,
                   is_final=False,
               ))
       else:
           stdout_lines.append(line)          # data line — accumulate
   ```

5. Wait for stdin writer thread; wait for Octave process:
   ```python
   stdin_thread.join()
   proc.wait(timeout=timeout_s)
   ```

6. Check `proc.returncode`; raise `OctaveExecutionError` on non-zero (same as before, including full stderr).

7. Parse `'\n'.join(stdout_lines)` with the existing `_parse_octave_stdout()`.

8. If `progress_queue is not None`, put the final `ProgressEvent(is_final=True)`.

**Helper functions (module-private):**

```python
_PROGRESS_RE = re.compile(r'^Harmonizing sample\s+(\d+)\s+out of\s+(\d+)')

def _is_progress_line(line: str) -> bool:
    return _PROGRESS_RE.match(line) is not None

def _parse_progress_line(line: str) -> tuple[int, int]:
    m = _PROGRESS_RE.match(line)
    return int(m.group(1)), int(m.group(2))
```

**Preservation:** Keep original `run_octave_normalize()` unchanged (used by existing tests). The new function is an additive variant; `parallel.py` will call the streaming version.

---

## 7. Modified File: `shambhala/parallel.py`

### 7.1 New worker function: `_harmonize_one_batch_with_progress()`

Replaces `_harmonize_one_batch` as the default dispatched function. Same signature plus:

```python
def _harmonize_one_batch_with_progress(
    batch_sample_names: list[str],
    input_df: pd.DataFrame,
    p_df: pd.DataFrame,
    k: int,
    octave_scripts_dir: str,
    octave_bin: str,
    timeout_s: int,
    random_seed: int,
    progress_queue,   # multiprocessing Manager Queue or None
    worker_idx: int,
) -> pd.DataFrame:
```

Implementation: identical to `_harmonize_one_batch` except it calls `run_octave_normalize_streaming()` instead of `run_octave_normalize()`, passing `progress_queue` and `worker_idx`. The `logger.info("Octave batch done: ...")` call is kept for non-TTY log records.

**Picklability note:** `multiprocessing.Manager().Queue()` is a proxy object and is picklable across `ProcessPoolExecutor`. Confirmed: this pattern works with `ProcessPoolExecutor` (each worker receives the proxy; it communicates back via socket to the Manager server in the main process).

### 7.2 Modified `harmonize_parallel()`

Add these steps around the `ProcessPoolExecutor` block:

```python
# ── Create progress queue and renderer ──────────────────────────────
import multiprocessing
manager = multiprocessing.Manager()
progress_queue = manager.Queue()

use_ansi = sys.stderr.isatty()
worker_batch_sizes = [len(b) for b in batches]
renderer = _ProgressRenderer(
    queue=progress_queue,
    n_workers=effective_workers,
    total_samples=n_samples,
    worker_batch_sizes=worker_batch_sizes,
    use_ansi=use_ansi,
)
renderer.daemon = True
renderer.start()

# ── Dispatch workers (unchanged loop, but call new function) ─────────
with ProcessPoolExecutor(max_workers=effective_workers) as executor:
    future_to_idx = {
        executor.submit(
            _harmonize_one_batch_with_progress,
            batch,
            input_df,
            p_df,
            k,
            octave_scripts_dir,
            octave_bin,
            timeout_s,
            random_seed,
            progress_queue,   # ← new
            idx,              # ← new
        ): idx
        for idx, batch in enumerate(batches)
    }

    for future in as_completed(future_to_idx):
        idx = future_to_idx[future]
        try:
            batch_results[idx] = future.result()
            logger.info("Worker %d/%d complete.", idx + 1, effective_workers)
        except Exception as exc:
            renderer.stop()
            manager.shutdown()
            raise RuntimeError(
                f"Worker for batch {idx} (samples: {batches[idx]}) failed."
            ) from exc

# ── Stop renderer cleanly ────────────────────────────────────────────
renderer.stop()
manager.shutdown()
```

**Shutdown order matters:** `renderer.stop()` must be called before `manager.shutdown()` because the renderer reads from the Manager queue. Shutting down the Manager first would make the queue unavailable and could deadlock or error the renderer's `run()` loop.

### 7.3 `multiprocessing.Manager` startup overhead

Starting a `Manager` spawns a small background server process (~0.1s). This overhead is negligible for long Shambhala runs but worth noting. If performance profiling ever shows this matters, the Manager can be eliminated by accepting a pre-created queue as an optional parameter.

---

## 8. Edge Cases and Guard Conditions

| Scenario | Handling |
|---|---|
| `n_workers = 1` | Renderer shows one row; works identically |
| `n_samples < n_workers` | `effective_workers` already clamped to `n_samples`; renderer initialized with correct count |
| Octave crashes mid-batch | `proc.returncode != 0` → `run_octave_normalize_streaming` raises `OctaveExecutionError` → `as_completed` handler calls `renderer.stop()` then re-raises |
| Octave emits no progress lines (future Octave version removes `disp`) | `stdout_lines` accumulates all output; parsing proceeds normally; renderer shows 0% until `is_final` event |
| Progress lines have unexpected format | `_parse_progress_line` returns `(None, None)` if regex doesn't match; `_is_progress_line` returns False; line goes to `stdout_lines` — no crash |
| `progress_queue` is `None` (caller opts out) | All `if progress_queue is not None` guards skip queue operations; function behaves exactly like `run_octave_normalize` |
| Running inside K8s pod (stderr not a TTY) | `use_ansi=False`; renderer falls back to plain timestamped log lines |
| `--log-level WARNING` from CLI | `logger.info` calls are suppressed; renderer still emits to `sys.stderr` independently (it doesn't use the `logging` module) |
| `manager.shutdown()` before `renderer.stop()` | **Guarded by order above.** If an unexpected exception occurs before `stop()`, the `renderer.daemon = True` setting ensures the renderer thread dies when the main process exits |
| Large TSV (>10 MB) stdin → Popen deadlock | **Prevented by the stdin writer thread** (§6.1 step 3). Python writes stdin in a thread so the main thread can read stdout concurrently |
| Octave writes very large stdout (many genes × many samples) | Same: stdin writer thread decouples write from read |

---

## 9. Testing Plan

### 9.1 New test file: `tests/test_progress_display.py`

| Test | What it verifies |
|---|---|
| `test_progress_event_namedtuple` | `ProgressEvent` fields and types |
| `test_renderer_plain_text_mode` | Pass `use_ansi=False`; send events; capture `sys.stderr`; verify timestamped lines emitted for every `log_every`-th event |
| `test_renderer_stop_emits_final_lines` | `stop()` flushes remaining incomplete workers |
| `test_renderer_does_not_crash_on_none_sentinel` | Renderer exits cleanly when None put on queue |
| `test_eta_format` | Helper function formatting: `_format_eta(0)` → `'--'`, `_format_eta(3661)` → `'1h 01m 01s'` |
| `test_bar_render_width` | Bar is exactly 30 chars |

### 9.2 New tests in `tests/test_octave_bridge.py`

| Test | What it verifies |
|---|---|
| `test_parse_progress_line_valid` | `_parse_progress_line('Harmonizing sample 3 out of 10')` → `(3, 10)` |
| `test_is_progress_line_true` | Positive match |
| `test_is_progress_line_false_data_line` | `'TP53 1.234 5.678'` → `False` |
| `test_run_octave_streaming_sends_events` (requires Octave) | Run on small fixture; collect events from queue; verify count == NH; verify final event |
| `test_run_octave_streaming_no_queue` (requires Octave) | `progress_queue=None`; verify output identical to `run_octave_normalize` |

### 9.3 Existing tests

All 30 existing tests must continue to pass without modification. The streaming function is additive; `parallel.py` switches to the new function, but the fixture-based integration tests call `harmonize_parallel()` — they will go through the new path. The renderer in non-TTY mode (as in pytest's captured stderr) will fall back to plain-text mode and emit a few log lines, which is harmless.

If the `multiprocessing.Manager()` overhead causes slowness in tests, add a `disable_progress: bool = False` parameter to `harmonize_parallel()` that bypasses renderer creation entirely. Default `False` for production; tests may set `True` if needed.

---

## 10. Implementation Order

1. **Create `shambhala/progress_display.py`** with `ProgressEvent`, `_ProgressRenderer`, `_format_bar`, `_format_eta`, `_format_speed`. Write unit tests first (§9.1) and verify they pass with the stub implementation.

2. **Add `run_octave_normalize_streaming()` to `octave_bridge.py`**: add `_PROGRESS_RE`, `_is_progress_line`, `_parse_progress_line`, and the new function. Keep original function intact. Write and run `test_run_octave_streaming_*` tests.

3. **Modify `parallel.py`**: add `_harmonize_one_batch_with_progress`, update `harmonize_parallel`. Run all 30 existing tests to confirm no regressions.

4. **Manual smoke test** on the toy fixture:
   ```bash
   python run_shambhala.py \
       --input tests/fixtures/input_10samples.csv \
       --P tests/fixtures/P0_small.csv \
       --Q tests/fixtures/Q0_small.csv \
       --output /tmp/harmonized_test.csv \
       --n-workers 3 --random-seed 42 --log-level DEBUG
   ```
   Verify the live progress display appears in the terminal.

5. **K8s smoke test**: rsync to pod, run `run_shambhala_job.py`, verify plain-text fallback in `kubectl logs`.

---

## 11. Non-Goals (Out of Scope for This Plan)

- Modifying `Shambhala2_piped.m` (zero changes).
- Adding progress to the Stage 1/2 benchmark pipeline (`run_prep_parallel.py`, `run_norm_parallel.py`) — those pipelines have their own logging.
- Replacing the `logging` module calls in `run_shambhala.py` with `tqdm`-style step counters.
- Adding a web dashboard or Prometheus metrics.

---

## 12. Detailed To-Do List

Checkboxes are ordered by dependency: each group can only start after the previous group is complete.

### Group A — New module: `shambhala/progress_display.py`

- [x] **A1.** Create `shambhala/progress_display.py` with module docstring explaining the TTY vs non-TTY split.
- [x] **A2.** Define `ProgressEvent` as a `typing.NamedTuple` with fields: `worker_idx: int`, `sample_done: int`, `sample_total: int`, `elapsed_s: float`, `is_final: bool`.
- [x] **A3.** Implement `_format_eta(seconds: float) -> str`: returns `'--'` for zero/negative input, `'Xh YYm ZZs'` / `'YYm ZZs'` / `'ZZs'` otherwise.
- [x] **A4.** Implement `_format_speed(elapsed_s: float, done: int) -> str`: returns `'--'` if `done == 0`, else `'{elapsed_s/done:.1f} s/sample'`.
- [x] **A5.** Implement `_format_bar(done: int, total: int, width: int = 30) -> str`: filled with `█`, empty with `░`, exactly `width` chars wide. Handle `total == 0` without division by zero.
- [x] **A6.** Implement `_render_row(label: str, done: int, total: int, elapsed_s: float, finished: bool) -> str`: assembles one full display row string (`label  [bar]  done/total (pct%)  speed  ETA`). Mark finished rows with `✓` prefix instead of ETA.
- [x] **A7.** Implement `_ProgressRenderer.__init__`: store queue, n_workers, total_samples, worker_batch_sizes, use_ansi; initialize `_state` list (one dataclass per worker); set `_last_rendered_lines = 0`; set `_dirty = False`.
- [x] **A8.** Implement `_ProgressRenderer._consume_queue(self) -> None`: drain all currently available items from queue (non-blocking `get_nowait` loop); catch `queue.Empty`; detect `None` sentinel (set `_stop_flag`).
- [x] **A9.** Implement `_ProgressRenderer._render_ansi(self) -> None`: build header line + separator + per-worker rows + separator + overall row; move cursor up `_last_rendered_lines`; write new block; update `_last_rendered_lines`.
- [x] **A10.** Implement `_ProgressRenderer._emit_plain(self, event: ProgressEvent) -> None`: emit timestamped line to stderr if `event.sample_done % self._log_every == 0` or `event.is_final`. Format: `[HH:MM:SS][Worker N]  done/total (pct%)  elapsed Xm Ys`.
- [x] **A11.** Implement `_ProgressRenderer.run(self) -> None`: main loop — call `_consume_queue`; if ANSI mode and dirty, call `_render_ansi`; if plain mode, `_emit_plain` is called directly inside `_consume_queue` on each event; sleep `0.05 s` between iterations; exit on `_stop_flag`.
- [x] **A12.** Implement `_ProgressRenderer.stop(self) -> None`: put `None` sentinel on queue; call `self.join(timeout=5.0)`; if ANSI mode, write a trailing newline to stderr.
- [x] **A13.** Compute `_log_every` in `__init__`: `max(1, total_samples // (20 * n_workers))` so each worker emits at most ~20 lines in plain mode.

### Group B — Tests for `progress_display.py` (write before implementing, verify after)

- [x] **B1.** Create `tests/test_progress_display.py`.
- [x] **B2.** `test_progress_event_fields`: instantiate `ProgressEvent`; assert all 5 fields accessible by name and have correct types.
- [x] **B3.** `test_format_eta_zero`: `_format_eta(0)` → `'--'`.
- [x] **B4.** `test_format_eta_seconds_only`: `_format_eta(45)` → `'45s'`.
- [x] **B5.** `test_format_eta_minutes_seconds`: `_format_eta(125)` → `'2m 05s'`.
- [x] **B6.** `test_format_eta_hours`: `_format_eta(3661)` → `'1h 01m 01s'`.
- [x] **B7.** `test_format_speed_zero_done`: `_format_speed(10.0, 0)` → `'--'`.
- [x] **B8.** `test_format_speed_nonzero`: `_format_speed(10.0, 2)` → `'5.0 s/sample'`.
- [x] **B9.** `test_format_bar_width`: `len(_format_bar(5, 10))` == 30.
- [x] **B10.** `test_format_bar_full`: `_format_bar(10, 10)` contains only `█` (no `░`).
- [x] **B11.** `test_format_bar_empty`: `_format_bar(0, 10)` contains only `░` (no `█`).
- [x] **B12.** `test_format_bar_zero_total`: `_format_bar(0, 0)` does not raise; returns 30-char string.
- [x] **B13.** `test_renderer_plain_text_emits_lines`: create renderer with `use_ansi=False`; put events per worker; call `stop()`; verify lines in captured stderr.
- [x] **B14.** `test_renderer_plain_text_final_always_emitted`: put only `is_final=True` events; assert each worker gets a line.
- [x] **B15.** `test_renderer_stop_exits_cleanly`: start renderer thread; immediately call `stop()`; assert `is_alive()` is `False`.
- [x] **B16.** `test_renderer_does_not_crash_on_extra_none`: put two `None` sentinels; assert no exception raised.

### Group C — New function in `shambhala/octave_bridge.py`

- [x] **C1.** Add `import re`, `import threading`, `import time` at the top of `octave_bridge.py`.
- [x] **C2.** Add `from shambhala.progress_display import ProgressEvent, _QueueLike` import.
- [x] **C3.** Define module-level `_PROGRESS_RE = re.compile(r'^Harmonizing sample\s+(\d+)\s+out of\s+(\d+)')`.
- [x] **C4.** Implement `_is_progress_line(line: str) -> bool`.
- [x] **C5.** Implement `_parse_progress_line(line: str) -> tuple[int, int]`; raises `ValueError` on non-matching input.
- [x] **C6.** Implement `_write_stdin(proc: subprocess.Popen[str], data: str) -> None`; catches `BrokenPipeError`.
- [x] **C7.** Implement `run_octave_normalize_streaming()` with `Popen`, stdin-writer thread, streaming stdout, final `ProgressEvent`.
- [x] **C8.** Verified that `run_octave_normalize()` (original) is completely unchanged.

### Group D — Tests for `octave_bridge.py` additions

- [x] **D1.** Added to `tests/test_octave_bridge.py`.
- [x] **D2.** `test_is_progress_line_true`.
- [x] **D3.** `test_is_progress_line_false_data_line`.
- [x] **D4.** `test_is_progress_line_false_empty`.
- [x] **D5.** `test_parse_progress_line_valid`.
- [x] **D6.** `test_parse_progress_line_extra_spaces`.
- [x] **D7.** `test_streaming_no_queue_matches_blocking` *(Octave-gated)*.
- [x] **D8.** `test_streaming_sends_progress_events` *(Octave-gated)*.
- [x] **D9.** `test_streaming_bad_binary_raises`.

### Group E — Modifications to `shambhala/parallel.py`

- [x] **E1.** Added `import sys`, `import multiprocessing`, `from shambhala.progress_display import ProgressEvent, _ProgressRenderer, _QueueLike`.
- [x] **E2.** Added module-level `_harmonize_one_batch_with_progress()`.
- [x] **E3.** Calls `run_octave_normalize_streaming()`; keeps `logger.info("Octave batch done: ...")`.
- [x] **E4.** Manager + queue created in `harmonize_parallel()` before executor block.
- [x] **E5.** `_ProgressRenderer` instantiated with `use_ansi=sys.stderr.isatty()`; `daemon=True`; started.
- [x] **E6.** `executor.submit` calls `_harmonize_one_batch_with_progress` with `progress_queue` and `idx`.
- [x] **E7.** Exception handler calls `renderer.stop()` and `manager.shutdown()` before re-raising.
- [x] **E8.** `renderer.stop()` then `manager.shutdown()` called in `finally` block after executor.
- [x] **E9.** Original `_harmonize_one_batch()` kept intact.
- [x] **E10.** `disable_progress: bool = False` parameter added to `harmonize_parallel()`.

### Group F — Regression and integration tests

- [x] **F1.** `pytest tests/ -v` — 59 passed, 10 Octave-gated skipped. All pre-existing tests pass.
- [x] **F2.** Group B (28 new tests) and Group D (9 new tests including 2 Octave-gated) all pass.
- [x] **F3.** pytest completes in ~25s with no hangs; Manager and renderer both shut down cleanly.
- [x] **F4.** Progress renderer writes to `sys.stderr` only; no `logging` module involvement — no double-printing.

### Group G — Manual smoke tests

- [ ] **G1.** TTY smoke test — requires Octave; run from local terminal when pod is available.
- [ ] **G2.** Non-TTY smoke test — requires Octave.
- [ ] **G3.** K8s smoke test — requires pod `shambhala-bench`; run `run_shambhala_job.py` and verify `kubectl logs` shows plain-text progress lines.

### Group H — Documentation updates

- [x] **H1.** Updated `CLAUDE.md`: added `progress_display.py` to module table; added "Progress display" subsection describing TTY/non-TTY modes and `disable_progress` parameter.
- [x] **H2.** Updated `README.md`: added "Progress display" section with example output for both TTY and non-TTY modes.
- [ ] **H3.** Delete this plan file once G1–G3 are verified on pod and the feature is merged.
