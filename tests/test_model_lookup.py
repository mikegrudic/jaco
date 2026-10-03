"""jaco.models.<name> is the model's subpackage whatever was imported first."""

import os
import subprocess
import sys
from pathlib import Path

import jaco

MODELS = ("starforge", "starforge_legacy", "wind_comparison")


def run_fresh(code):
    """Run code in a fresh interpreter importing this jaco"""
    env = dict(os.environ)
    env["PYTHONPATH"] = os.pathsep.join([str(Path(jaco.__file__).resolve().parents[1])]
                                        + [p for p in [env.get("PYTHONPATH")] if p])
    return subprocess.run([sys.executable, "-W", "ignore", "-c", code], env=env, capture_output=True, text=True)


def test_from_models_import_gives_subpackages():
    check = "; ".join(f"assert isinstance({m}, types.ModuleType) and hasattr({m}, 'make_model'), {m}" for m in MODELS)
    result = run_fresh(f"import types\nfrom jaco.models import {', '.join(MODELS)}\n{check}")
    assert result.returncode == 0, result.stderr


def test_attribute_lookup_does_not_depend_on_import_order():
    result = run_fresh("import importlib, jaco.models\n"
                       "first = getattr(jaco.models, 'starforge', None)\n"
                       "module = importlib.import_module('jaco.models.starforge')\n"
                       "assert first is None or first is module, first\n"
                       "assert jaco.models.starforge is module")
    assert result.returncode == 0, result.stderr
