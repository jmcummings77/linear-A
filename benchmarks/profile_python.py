#!/usr/bin/env python3
"""Sample only the main thread of the Python runner executed in this process."""
import json
from pathlib import Path
import runpy
import sys
import threading
import time


def main():
    if len(sys.argv) < 3:
        raise SystemExit("usage: profile_python.py OUTPUT RUNNER [ARGUMENTS...]")
    output = Path(sys.argv[1])
    runner = Path(sys.argv[2]).resolve()
    arguments = sys.argv[3:]
    samples = []
    done = threading.Event()
    main_thread = threading.get_ident()

    def sample():
        while not done.wait(0.001):
            frame = sys._current_frames().get(main_thread)
            stack = []
            while frame is not None:
                stack.append("%s (%s:%s)" % (frame.f_code.co_name, Path(frame.f_code.co_filename).name, frame.f_lineno))
                frame = frame.f_back
            if stack:
                samples.append({"at_ns": time.perf_counter_ns(), "stack": stack[::-1]})

    worker = threading.Thread(target=sample, name="matrix-stack-sampler", daemon=True)
    sys.path.insert(0, str(runner.parent))
    sys.argv = [str(runner)] + arguments
    worker.start()
    try:
        runpy.run_path(str(runner), run_name="__main__")
    finally:
        end_ns = time.perf_counter_ns()
        done.set()
        worker.join(timeout=2)
        output.write_text(json.dumps({"samples": samples, "end_ns": end_ns}), encoding="utf-8")


if __name__ == "__main__":
    main()
