API Documentation
=================

Process
-------

The ``Process`` class is the top-level building block for assembling microphysics
models. Processes represent individual physical mechanisms (ionization,
recombination, cooling, heating) and can be composed with ``+`` to build
composite networks.

.. autoclass:: jaco.Process
   :members:
   :undoc-members:

Processes are immutable. Physics is written as reactions and thermal terms; the other names are thin factories
returning them:

.. automodule:: jaco.processes
   :members: Reaction, ThermalTerm, collisional_thermal_term, Ionization, Recombination, CollisionalIonization, GasPhaseRecombination, FreeFreeEmission, LineCoolingSimple
   :undoc-members:

Model
-----

A ``Model`` is a keyed collection of processes (``+`` adds one and raises on a duplicate id; ``without(id)``,
``replace(id, new)``) plus the declarations that make their rate equations solvable: solve variables, the
time-dependent set, steady-state species (closed from their own rate equation when the network is assembled), fixed
abundances, derived parameters, intermediates and model-level ``Rule`` rewrites of the processes. Everything is
worked out from the current processes when the network is assembled, so ``model + process`` is the model with that
process in every respect.

.. autoclass:: jaco.model.Model
   :members:

.. autoclass:: jaco.model.Rule

Declarations: every species with its kind, and every parameter besides the core ones and those the species imply.
Code generation refuses undeclared symbols, and the header defines ``JACO_HAS_VAR_<name>``/``JACO_HAS_PARAM_<name>``.

.. automodule:: jaco.declarations
   :members: Parameter, Species, CORE_PARAMETERS

EquationSystem
--------------

The ``EquationSystem`` is the core symbolic engine that stores the rate
equations for all species and provides methods for reducing, solving, and
generating code from the network. A model's ``network`` is one, with the
model's declarations attached.

.. autoclass:: jaco.EquationSystem
   :members:
   :undoc-members:

- **Reduction pipeline** (``reduced()``):

  1. Substitute derived parameters
  2. Set time dependence (BDF for evolved species, steady-state for others)
  3. Conservation reductions (n->x conversion, charge neutrality, atom conservation)
  4. Fix species (substitute fixed and steady-state abundances, remove their equations)
  5. Prune decoupled equations (remove orphaned equations like dust heat)

  Every removed equation is recorded with its reason in ``discarded``; code generation prints them.

- **Code generation** (``generate_code()``): Produces C/C++/CUDA/Python/Julia
  source files with the RHS function, Jacobian, EOS functions, and
  interpolation table infrastructure.

Models
------

Each subpackage of ``jaco.models`` provides ``make_model()``, returning a ``Model``.

Starforge
^^^^^^^^^

The STARFORGE model implements ISM thermochemistry with H2, metals, dust,
cosmic rays, and radiation. It is the primary model for star formation
simulations in GIZMO. ``starforge_legacy`` builds GIZMO's legacy cooling
module from the same process library.

.. autofunction:: jaco.models.starforge.starforge.make_model

.. autofunction:: jaco.models.starforge_legacy.make_model

Wind Comparison
^^^^^^^^^^^^^^^

A simple atomic-H cooling model for wind bubble tests.

.. autofunction:: jaco.models.wind_comparison.cooling.make_model

GIZMO Code Generation
----------------------

The GIZMO-specific codegen layer writes source files for integration into GIZMO's
build system.

.. automodule:: jaco.codegen.gizmo
   :members: generate_funcjac_code

Model-specific defaults for solve variables and time-dependent species:

.. data:: jaco.codegen.gizmo._MODEL_DEFAULTS

   Dictionary mapping model names to their default ``solve_vars`` and
   ``time_dependent`` lists. Used when the codegen is invoked from the
   command line (``python -m jaco.codegen.gizmo <model_name>``).

Equation
--------

.. autoclass:: jaco.Equation
   :members:
   :undoc-members:

Numerical Solvers
-----------------

.. automodule:: jaco.numerics
   :members: newton_rootsolve

EOS
---

.. automodule:: jaco.eos.eos
   :members:
   :undoc-members:

Symbols
-------

.. automodule:: jaco.symbols
   :members: T, n_, x_, n_Htot, sanitize_name, sanitize_symbols, table_interp_2d
   :undoc-members:

Interpolation
-------------

.. automodule:: jaco.interpolation
   :members: PiecewiseLinearInterp, PiecewiseConstantInterp, TableInterp2D, register_table, save_table_hdf5
   :undoc-members:
