"""Implementation of the Model class: comprises a set of equations, symbols, and default assumptions that together
constitute a microphysics model.
"""

from dataclasses import dataclass
from ..process import Process
from ..equation_system import EquationSystem


@dataclass
class Model:
    processes: list  # list of processes
    assumptions: dict  # mappings for variable substitutions when not specified - can be symbolic or numerical - specify the EOS at this level?

    @property
    def process(self) -> Process:
        """Returns the full, summed list of processes constituting the system of equations."""
        return sum(self.processes)

    @property
    def network(self) -> EquationSystem:
        """Returns the system of equations for the processes."""
        return self.process.network

    @property
    def bibliography(self) -> dict:
        """Returns a dict of processes and respective bibcodes."""
        return {p.name: p.bibliography for p in self.processes}
