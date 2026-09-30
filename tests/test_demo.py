"""Smoke test: the README quickstart demo runs and prints what the README shows."""

import runpy
from pathlib import Path

DEMO = Path(__file__).resolve().parent.parent / "examples" / "demo.py"


def test_demo_output(capsys):
    runpy.run_path(str(DEMO), run_name="__main__")
    out = capsys.readouterr().out
    assert "SNPs on at least one array: 3 of 4" in out
    assert "rs9100001  ARRAY_A,ARRAY_B   10000   225  2550     2775" in out
    assert "rs9100009  -                     0     0     0        0" in out
