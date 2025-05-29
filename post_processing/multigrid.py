import numpy as np

def laplacian(f, dx):
    return (
        np.roll(f, 1, axis=0) + np.roll(f, -1, axis=0) +
        np.roll(f, 1, axis=1) + np.roll(f, -1, axis=1) - 4*f
    ) / dx**2

def smooth(psi, rhs, dx, iterations=3):
    """
    Using Gauss-Seidel method as pre- and post-smoothing to reduce high frequency
    errors.
    """
    for _ in range(iterations):
        psi = 0.25 * (
            np.roll(psi, 1, axis=0) + np.roll(psi, -1, axis=0) +
            np.roll(psi, 1, axis=1) + np.roll(psi, -1, axis=1) -
            dx**2 * rhs
        )
    return psi

def residual(psi, rhs, dx):
    # Compute residual errors ：res = rhs - Lψ
    return rhs - laplacian(psi, dx)

def restrict(fine):
    """
    Downsampling (restricting) the residual error to a coarser grid.
    """
    Nc = fine.shape[0] // 2
    coarse = np.zeros((Nc, Nc))
    for i in range(Nc):
        for j in range(Nc):
            ii, jj = 2 * i, 2 * j
            coarse[i, j] = 0.25 * (fine[ii, jj] + fine[ii+1, jj] +
                                   fine[ii, jj+1] + fine[ii+1, jj+1])
    return coarse

def prolong(coarse):
    """
    Interpolating (prolonging) a correction computed on a coarser grid into a finer grid.
    """
    Nc = coarse.shape[0]
    Nf = Nc * 2
    fine = np.zeros((Nf, Nf))
    for i in range(Nc):
        for j in range(Nc):
            fine[2*i:2*i+2, 2*j:2*j+2] = coarse[i, j]
    return fine

def v_cycle(rhs, dx, level):
    """
    V-cycle Multigrid main function
    - level = 0 means the coarsest gird
    - Otherwise, recursively downsampling, computing errors and prolonging
    """
    N = rhs.shape[0]
    psi = np.zeros_like(rhs)

    # Stop recursion at smallest grid size, otherwise continue recursion
    if level == 0 or N <= 4:
        return smooth(psi, rhs, dx, iterations=50)

    # Pre-Smoothing
    psi = smooth(psi, rhs, dx)
    # Compute Residual Errors
    res = residual(psi, rhs, dx)
    # Restriction
    res_c = restrict(res)
    # Recursively compute errors
    err_c = v_cycle(res_c, 2*dx, level-1)
    # Prolongation and Correction
    err_f = prolong(err_c)
    psi += err_f
    # Post-Smoothing
    psi = smooth(psi, rhs, dx)

    return psi

def w_cycle(rhs, dx, level):
    """
    W-cycle Multigrid main function
    - level = 0 means the coarsest grid
    - Otherwise, recursively apply two coarse-level corrections
    """
    N = rhs.shape[0]
    psi = np.zeros_like(rhs)

    # Stop recursion at coarsest grid
    if level == 0 or N <= 4:
        return smooth(psi, rhs, dx, iterations=50)

    # Pre-smoothing
    psi = smooth(psi, rhs, dx)

    # Compute residual
    res = residual(psi, rhs, dx)

    # Restrict to coarse grid
    res_c = restrict(res)

    # First coarse-grid correction
    err_c1 = w_cycle(res_c, 2*dx, level - 1)

    # Second coarse-grid correction (W-cycle key step)
    err_c2 = w_cycle(err_c1, 2*dx, level - 1)

    # Prolongate and apply correction
    err_f = prolong(err_c2)
    psi += err_f

    # Post-smoothing
    psi = smooth(psi, rhs, dx)

    return psi

def solve_poisson(zeta, dx):
    # Given ζ，solve Δψ = ζ
    rhs = zeta
    levels = int(np.log2(zeta.shape[0])) - 2
    return v_cycle(rhs, dx, levels)
    #return w_cycle(rhs, dx, levels)