# from jaco.models.starforge.gas_dust_collisions import GasDustCollisions
from ..gas_dust_collisions import GasDustCollisions


def test_dust_heat():
    """Simple check that the process conserves energy"""
    gas_dust_collisions = GasDustCollisions()
    assert gas_dust_collisions.network["dust heat"].rhs == -gas_dust_collisions.network["heat"].rhs
    assert gas_dust_collisions.heat == -gas_dust_collisions.dust_heat
