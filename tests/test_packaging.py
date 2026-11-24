"""The runtime dependencies in pyproject.toml are exactly the third-party packages jaco imports."""

import ast
import re
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]


def test_runtime_dependencies_match_imports():
    tomllib = pytest.importorskip("tomllib")
    deps = tomllib.loads((ROOT / "pyproject.toml").read_text())["project"]["dependencies"]
    declared = {re.split(r"[\s<>=!~\[;]", d, maxsplit=1)[0].lower() for d in deps}

    package = ROOT / "src" / "jaco"
    imported = set()
    for path in package.rglob("*.py"):
        if "tests" in path.relative_to(package).parts or path.name.startswith("test_"):
            continue
        for node in ast.walk(ast.parse(path.read_text())):
            if isinstance(node, ast.Import):
                imported |= {a.name.split(".")[0] for a in node.names}
            elif isinstance(node, ast.ImportFrom) and node.level == 0:
                imported.add(node.module.split(".")[0])
    third_party = imported - set(sys.stdlib_module_names) - {"jaco"}
    assert declared == third_party, (f"imported but not declared: {sorted(third_party - declared)}; "
                                     f"declared but not imported: {sorted(declared - third_party)}")
