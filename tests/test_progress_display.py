"""
Unit tests for shambhala/progress_display.py.
"""

import io
import pathlib
import queue
import sys
import time

import pytest

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent.parent))

from shambhala.progress_display import (
    ProgressEvent,
    _ProgressRenderer,
    _format_bar,
    _format_eta,
    _format_speed,
    _render_row,
)


def _make_renderer(
    n_workers: int = 2,
    total_samples: int = 20,
    batch_sizes: list[int] | None = None,
    use_ansi: bool = False,
) -> tuple[_ProgressRenderer, queue.Queue[ProgressEvent | None], io.StringIO]:
    q: queue.Queue[ProgressEvent | None] = queue.Queue()
    fake_stderr = io.StringIO()
    if batch_sizes is None:
        batch_sizes = [total_samples // n_workers] * n_workers
    renderer = _ProgressRenderer(
        queue=q,
        n_workers=n_workers,
        total_samples=total_samples,
        worker_batch_sizes=batch_sizes,
        use_ansi=use_ansi,
        _stderr=fake_stderr,
    )
    renderer.daemon = True
    return renderer, q, fake_stderr


# ── ProgressEvent ──────────────────────────────────────────────────────────────


def test_progress_event_fields():
    ev = ProgressEvent(worker_idx=1, sample_done=5, sample_total=10, elapsed_s=2.5, is_final=False)
    assert ev.worker_idx == 1
    assert ev.sample_done == 5
    assert ev.sample_total == 10
    assert ev.elapsed_s == 2.5
    assert ev.is_final is False


def test_progress_event_is_namedtuple():
    ev = ProgressEvent(0, 3, 10, 1.0, True)
    assert isinstance(ev, tuple)
    assert len(ev) == 5


# ── _format_eta ────────────────────────────────────────────────────────────────


def test_format_eta_zero():
    assert _format_eta(0) == "--"


def test_format_eta_negative():
    assert _format_eta(-5.0) == "--"


def test_format_eta_seconds_only():
    assert _format_eta(45) == "45s"


def test_format_eta_minutes_seconds():
    assert _format_eta(125) == "2m 05s"


def test_format_eta_hours():
    assert _format_eta(3661) == "1h 01m 01s"


def test_format_eta_exact_minute():
    assert _format_eta(60) == "1m 00s"


def test_format_eta_exact_hour():
    assert _format_eta(3600) == "1h 00m 00s"


# ── _format_speed ──────────────────────────────────────────────────────────────


def test_format_speed_zero_done():
    assert _format_speed(10.0, 0) == "--"


def test_format_speed_nonzero():
    assert _format_speed(10.0, 2) == "5.0 s/sample"


def test_format_speed_one_sample():
    assert _format_speed(3.7, 1) == "3.7 s/sample"


# ── _format_bar ────────────────────────────────────────────────────────────────


def test_format_bar_width():
    assert len(_format_bar(5, 10)) == 30


def test_format_bar_custom_width():
    assert len(_format_bar(5, 10, width=20)) == 20


def test_format_bar_full():
    bar = _format_bar(10, 10)
    assert "░" not in bar
    assert len(bar) == 30


def test_format_bar_empty():
    bar = _format_bar(0, 10)
    assert "█" not in bar
    assert len(bar) == 30


def test_format_bar_zero_total():
    bar = _format_bar(0, 0)
    assert len(bar) == 30
    assert "█" not in bar


def test_format_bar_half():
    bar = _format_bar(5, 10)
    assert bar.count("█") == 15
    assert bar.count("░") == 15


# ── _render_row ────────────────────────────────────────────────────────────────


def test_render_row_contains_bar():
    row = _render_row("Worker 1", 5, 10, 10.0, False)
    assert "[" in row and "]" in row
    assert "5/10" in row


def test_render_row_finished_shows_done():
    row = _render_row("Worker 1", 10, 10, 20.0, True)
    assert "✓ done" in row
    assert "ETA" not in row


def test_render_row_unfinished_shows_eta():
    row = _render_row("Worker 1", 5, 10, 10.0, False)
    assert "ETA" in row


def test_render_row_bar_width():
    row = _render_row("W", 5, 10, 10.0, False)
    start = row.index("[") + 1
    end = row.index("]")
    bar_content = row[start:end]
    assert len(bar_content) == 30


# ── _ProgressRenderer plain-text mode ─────────────────────────────────────────


def test_renderer_plain_text_emits_lines():
    renderer, q, fake_stderr = _make_renderer(n_workers=2, total_samples=20)
    renderer.start()
    for i in range(10):
        q.put(ProgressEvent(0, i + 1, 10, float(i + 1), i == 9))
    for i in range(10):
        q.put(ProgressEvent(1, i + 1, 10, float(i + 1), i == 9))
    renderer.stop()
    output = fake_stderr.getvalue()
    assert "[Worker 1]" in output or "[Worker 2]" in output


def test_renderer_plain_text_final_always_emitted():
    renderer, q, fake_stderr = _make_renderer(n_workers=2, total_samples=200, batch_sizes=[100, 100])
    renderer.start()
    # Only send final events, log_every would normally suppress intermediate lines
    q.put(ProgressEvent(0, 100, 100, 10.0, True))
    q.put(ProgressEvent(1, 100, 100, 10.0, True))
    renderer.stop()
    output = fake_stderr.getvalue()
    assert "[Worker 1]" in output
    assert "[Worker 2]" in output


def test_renderer_stop_exits_cleanly():
    renderer, q, _ = _make_renderer()
    renderer.start()
    renderer.stop()
    assert not renderer.is_alive()


def test_renderer_stop_after_empty_queue_exits():
    renderer, q, _ = _make_renderer()
    renderer.start()
    time.sleep(0.1)
    renderer.stop()
    assert not renderer.is_alive()


def test_renderer_does_not_crash_on_extra_none():
    renderer, q, _ = _make_renderer()
    renderer.start()
    q.put(None)
    q.put(None)
    # Should not deadlock or crash — join with generous timeout
    renderer.join(timeout=3.0)
    assert not renderer.is_alive()


def test_renderer_processes_all_events_before_stop():
    renderer, q, fake_stderr = _make_renderer(n_workers=1, total_samples=5, batch_sizes=[5])
    renderer.start()
    for i in range(5):
        q.put(ProgressEvent(0, i + 1, 5, float(i + 1), i == 4))
    renderer.stop()
    output = fake_stderr.getvalue()
    # Final event (5/5) must be logged
    assert "5/5" in output
