"""Command-line entry point: generate GIZMO's JACO microphysics sources for a jaco model.

    python -m jaco.codegen.gizmo.gizmo starforge --language c --ext .cc --output-dir cooling

The model is ``jaco.models.<model>``; it must provide ``make_model()``, returning a :class:`jaco.model.Model` (which
declares its solve variables), and may set ``GIZMO_FAMILY``, the name of the host-code family macro
JACO_FAMILY_<FAMILY> for models that share GIZMO's code paths (default: the model name).
"""

import argparse
import sys
from importlib import import_module


def _split(s):
    return [v for v in s.split(",") if v] if s else None


def main(argv=None):
    parser = argparse.ArgumentParser(description="Generate jaco microphysics sources for GIZMO")
    parser.add_argument("model", help="model name, i.e. a module jaco.models.<model> (e.g. starforge)")
    parser.add_argument("--language", default="c", help="c, c++, cuda, python or julia (default: c)")
    parser.add_argument("--ext", default=None, help="source file extension, e.g. .cc")
    parser.add_argument("--output-dir", default=".", help="directory for the generated files")
    parser.add_argument("--jac-mode", default="symbolic", help="symbolic or autodiff")
    parser.add_argument("--solve-vars", default=None, help="comma-separated override of the model's solve variables")
    parser.add_argument("--time-dependent", default=None,
                        help="comma-separated override of the model's time-dependent variables")
    args = parser.parse_args(argv)

    model_mod = import_module(f"jaco.models.{args.model}")
    if not hasattr(model_mod, "make_model"):
        sys.exit(f"jaco.models.{args.model} has no make_model()")
    from . import generate_funcjac_code

    generate_funcjac_code(
        model_mod.make_model(),
        solve_vars=_split(args.solve_vars),
        time_dependent=_split(args.time_dependent),
        language=args.language,
        source_ext=args.ext,
        output_dir=args.output_dir,
        jac_mode=args.jac_mode,
        model_name=args.model,
        model_family=getattr(model_mod, "GIZMO_FAMILY", args.model),
    )


if __name__ == "__main__":
    main()
