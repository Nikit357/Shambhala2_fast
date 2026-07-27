"""
Live per-worker progress display for the Shambhala2 harmonization pipeline.

In TTY mode: renders fixed rows (one per worker) using ANSI escape codes,
refreshed up to 4 Hz. In non-TTY mode (log files, K8s): emits plain
timestamped lines to stderr.
"""

import queue as queue_module
import sys
import threading
import time
from dataclasses import dataclass
from typing import NamedTuple, Protocol, TextIO


class ProgressEvent(NamedTuple):
    worker_idx: int
    sample_done: int
    sample_total: int
    elapsed_s: float
    is_final: bool


class _QueueLike(Protocol):
    def put(self, item: ProgressEvent | None) -> None: ...
    def get_nowait(self) -> ProgressEvent | None: ...


@dataclass
class _WorkerState:
    done: int = 0
    total: int = 0
    elapsed_s: float = 0.0
    finished: bool = False


def _format_eta(seconds: float) -> str:
    if seconds <= 0:
        return "--"
    s = int(seconds)
    h, rem = divmod(s, 3600)
    m, sec = divmod(rem, 60)
    if h > 0:
        return f"{h}h {m:02d}m {sec:02d}s"
    if m > 0:
        return f"{m}m {sec:02d}s"
    return f"{sec}s"


def _format_speed(elapsed_s: float, done: int) -> str:
    if done == 0:
        return "--"
    return f"{elapsed_s / done:.1f} s/sample"


def _format_bar(done: int, total: int, width: int = 30) -> str:
    if total <= 0:
        return "░" * width
    filled = min(int(width * done / total), width)
    return "█" * filled + "░" * (width - filled)


def _render_row(label: str, done: int, total: int, elapsed_s: float, finished: bool) -> str:
    bar = _format_bar(done, total)
    pct = (done / total * 100) if total > 0 else 0.0
    speed = _format_speed(elapsed_s, done)
    if finished:
        return f" {label:<10}  [{bar}]  {done:4d}/{total} (100%)  {speed}  ✓ done"
    remaining = total - done
    rate = done / elapsed_s if elapsed_s > 0 and done > 0 else 0.0
    eta_s = remaining / rate if rate > 0 else 0.0
    eta_str = _format_eta(eta_s)
    return f" {label:<10}  [{bar}]  {done:4d}/{total} ({pct:3.0f}%)  {speed}  ETA {eta_str}"


class _ProgressRenderer(threading.Thread):
    def __init__(
        self,
        queue: _QueueLike,
        n_workers: int,
        total_samples: int,
        worker_batch_sizes: list[int],
        use_ansi: bool,
        _stderr: TextIO | None = None,
    ) -> None:
        super().__init__()
        self._queue = queue
        self._n_workers = n_workers
        self._total_samples = total_samples
        self._use_ansi = use_ansi
        self._stderr: TextIO = _stderr if _stderr is not None else sys.stderr
        self._log_every = max(1, total_samples // (20 * n_workers))
        self._state: list[_WorkerState] = [
            _WorkerState(total=worker_batch_sizes[i]) for i in range(n_workers)
        ]
        self._last_rendered_lines: int = 0
        self._dirty: bool = False
        self._stop_flag: bool = False
        self._start_time: float = 0.0

    def _consume_queue(self) -> None:
        while True:
            try:
                item = self._queue.get_nowait()
            except queue_module.Empty:
                break
            if item is None:
                self._stop_flag = True
                break
            event: ProgressEvent = item
            w = event.worker_idx
            self._state[w].done = event.sample_done
            self._state[w].total = event.sample_total
            self._state[w].elapsed_s = event.elapsed_s
            if event.is_final:
                self._state[w].finished = True
            self._dirty = True
            if not self._use_ansi:
                self._emit_plain(event)

    def _emit_plain(self, event: ProgressEvent) -> None:
        if event.sample_done % self._log_every == 0 or event.is_final:
            ts = time.strftime("%H:%M:%S")
            label = f"Worker {event.worker_idx + 1}"
            pct = (event.sample_done / event.sample_total * 100) if event.sample_total > 0 else 0.0
            elapsed_str = _format_eta(event.elapsed_s)
            self._stderr.write(
                f"[{ts}][{label}]  {event.sample_done}/{event.sample_total}"
                f" ({pct:.0f}%)  elapsed {elapsed_str}\n"
            )
            self._stderr.flush()

    def _render_ansi(self) -> None:
        total_done = sum(s.done for s in self._state)
        wall_elapsed = time.time() - self._start_time if self._start_time > 0 else 0.0
        sep = "─" * 79
        lines: list[str] = [
            f"[Shambhala] {total_done}/{self._total_samples} samples | {self._n_workers} workers",
            sep,
        ]
        for i, s in enumerate(self._state):
            label = f"Worker {i + 1}"
            lines.append(_render_row(label, s.done, s.total, s.elapsed_s, s.finished))
        lines.append(sep)
        rate = total_done / wall_elapsed if wall_elapsed > 0 and total_done > 0 else 0.0
        remaining = self._total_samples - total_done
        eta_s = remaining / rate if rate > 0 else 0.0
        pct = (total_done / self._total_samples * 100) if self._total_samples > 0 else 0.0
        overall_bar = _format_bar(total_done, self._total_samples)
        lines.append(
            f" {'Overall':<10}  [{overall_bar}]  {total_done:4d}/{self._total_samples}"
            f" ({pct:3.0f}%)  elapsed {_format_eta(wall_elapsed)}  ETA {_format_eta(eta_s)}"
        )
        if self._last_rendered_lines > 0:
            self._stderr.write(f"\033[{self._last_rendered_lines}A\033[J")
        output = "\n".join(lines) + "\n"
        self._stderr.write(output)
        self._stderr.flush()
        self._last_rendered_lines = len(lines)
        self._dirty = False

    def run(self) -> None:
        self._start_time = time.time()
        _RENDER_INTERVAL = 0.25
        last_render = 0.0
        while True:
            self._consume_queue()
            now = time.time()
            if self._use_ansi and self._dirty and (now - last_render) >= _RENDER_INTERVAL:
                self._render_ansi()
                last_render = now
            if self._stop_flag:
                break
            time.sleep(0.05)
        # Final flush after stop signal
        self._consume_queue()
        if self._use_ansi and self._dirty:
            self._render_ansi()

    def stop(self) -> None:
        self._queue.put(None)
        self.join(timeout=5.0)
        if self._use_ansi:
            self._stderr.write("\n")
            self._stderr.flush()
