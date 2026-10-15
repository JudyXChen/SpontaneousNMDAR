"""Run this condition with the shared legacy runtime."""
from pathlib import Path
_runtime = Path(__file__).absolute().parents[3] / "shared/legacy/model.py"
exec(compile(_runtime.read_text(), str(__file__), "exec"))
