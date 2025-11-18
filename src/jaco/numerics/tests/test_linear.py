import numpy as np
import jax, jax.numpy as jnp
from ..solvers import newton_rootsolve


def test_newton_rootsolve_linear(N=10**5, dim=10):
    """Test: solve a linear system of dim variables with the Newton solver and check the solution.
    Solver should get the right answer to machine precision.
    """

    A = np.random.rand(N, dim, dim)
    cond = np.linalg.cond(A)
    A, cond = A[cond < 1e3], cond[cond < 1e3]  # just take the easy ones, we don't need to be fancy here
    N = len(A)
    b = np.random.rand(N, dim, 1)
    sol = np.linalg.solve(A, b)[:, :, 0]
    b = b[:, :, 0]

    @jax.jit
    def func(x, *params):
        """horrid function for mapping X, *params signature to linear residual Ax-b"""
        dim = x.shape[0]
        A = jnp.array(params)[: dim * dim].reshape((dim, dim))
        b = jnp.array(params[dim * dim : dim * (dim + 1)])
        return jnp.matmul(A, x) - b

    p = np.c_[A.reshape(N, dim * dim), b.reshape(N, dim)]
    guesses = np.copy(b)
    sol_newton = newton_rootsolve(func, guesses, p)  # [:,:,None]

    error = np.sum((sol_newton - sol) ** 2, axis=1) ** 0.5
    norm = np.sum(sol**2, axis=1)
    assert np.all(error < 1e-2 * norm)
