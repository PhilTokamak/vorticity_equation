import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

def laplacian(f, dx):
    return (
        np.roll(f, 1, axis=0) + np.roll(f, -1, axis=0) +
        np.roll(f, 1, axis=1) + np.roll(f, -1, axis=1) - 4*f
    ) / dx**2

def compute_velocity(psi, dx):
    u =   (np.roll(psi, -1, axis=1) - np.roll(psi, 1, axis=1)) / (2 * dx)  # ∂ψ/∂y
    v = - (np.roll(psi, -1, axis=0) - np.roll(psi, 1, axis=0)) / (2 * dx)  # -∂ψ/∂x
    return u, v

def arakawa_jacobian(z, psi, dx):
    # Arakawa Jacobian J_A = (jplusplus + jpluscros + jcrossplus) / 3
    jplusplus = (
        (np.roll(z, -1, axis=0) - np.roll(z, 1, axis=0)) *
        (np.roll(psi, -1, axis=1) - np.roll(psi, 1, axis=1)) -
        (np.roll(z, -1, axis=1) - np.roll(z, 1, axis=1)) *
        (np.roll(psi, -1, axis=0) - np.roll(psi, 1, axis=0))
    )

    jpluscross = (
        np.roll(z, -1, axis=0) * (
            np.roll(psi, -1, (0, 1)) - np.roll(psi, (-1, 1), (0, 1))
        ) -
        np.roll(z, 1, axis=0) * (
            np.roll(psi, (1, -1), (0, 1)) - np.roll(psi, 1, (0, 1))
        ) -
        np.roll(z, -1, axis=1) * (
            np.roll(psi, -1, (0, 1)) - np.roll(psi, (1, -1), (0, 1))
        ) +
        np.roll(z, 1, axis=1) * (
            np.roll(psi, (-1, 1), (0, 1)) - np.roll(psi, 1, (0, 1))
        )
    )

    jcrossplus = (
        np.roll(zeta, -1, (0, 1)) * (np.roll(psi, -1, 1) - np.roll(psi, -1, 0)) -
        np.roll(zeta, 1, (0, 1)) * (np.roll(psi, 1, 0) - np.roll(psi, 1, 1)) -
        np.roll(zeta, (1, -1), (0, 1)) * (np.roll(psi, -1, 1) - np.roll(psi, 1, 0)) +
        np.roll(zeta, (-1, 1), (0, 1)) * (np.roll(psi, -1, 0) - np.roll(psi, 1, 1))
    )

    return (jplusplus + jpluscross + jcrossplus) / (12 * dx**2)

# def solve_poisson(zeta, dx, tol=1e-4, max_iter=10000):
#     psi = np.zeros_like(zeta)
#     for _ in range(max_iter):
#         psi_new = 0.25 * (
#             np.roll(psi, 1, axis=0) + np.roll(psi, -1, axis=0) +
#             np.roll(psi, 1, axis=1) + np.roll(psi, -1, axis=1) +
#             dx**2 * zeta
#         )
#         if np.max(np.abs(psi_new - psi)) < tol:
#             break
#         psi = psi_new
#     return psi

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

def rhs(zeta, dx, nu):
    psi = solve_poisson(zeta, dx)
    jac = arakawa_jacobian(zeta, psi, dx)
    diff = nu * laplacian(zeta, dx)
    return jac + diff

def rk4_step(zeta, dt, dx, nu):
    k1 = rhs(zeta, dx, nu)
    k2 = rhs(zeta + 0.5 * dt * k1, dx, nu)
    k3 = rhs(zeta + 0.5 * dt * k2, dx, nu)
    k4 = rhs(zeta + dt * k3, dx, nu)
    return zeta + (dt/6)*(k1 + 2*k2 + 2*k3 + k4)

def initialize_vorticity(N, kind='vortex'):
    x = np.linspace(0, 2*np.pi, N, endpoint=False)
    y = np.linspace(0, 2*np.pi, N, endpoint=False)
    X, Y = np.meshgrid(x, y, indexing='ij')
    if kind == 'vortex':
        r2 = (X - np.pi)**2 + (Y - np.pi)**2
        # r3 = (X - 3 * np.pi / 2)**2 + (Y - np.pi)**2
        return np.exp(-r2 / 0.1) # - np.exp(-r3 / 0.1)
    elif kind == 'random':
        return np.random.randn(N, N) * 0.1
    else:
        raise ValueError("Unknown kind")


def create_animation(zeta_history, T, filename="vorticity.mp4", interval=50):
    fig, ax = plt.subplots()
    im = ax.imshow(zeta_history[0], origin='lower', cmap='seismic',
                   extent=[0, 2*np.pi, 0, 2*np.pi], animated=True)
    cbar = plt.colorbar(im, ax=ax)

    def update(frame):
        im.set_array(zeta_history[frame])
        ax.set_title("t = {:.3f}".format(T[frame]))
        return [im]

    ani = animation.FuncAnimation(fig, update, frames=len(zeta_history),
                                  interval=interval, blit=True)

    ani.save(filename, writer='ffmpeg', dpi=150)
    plt.close()

def create_velocity_animation(velocity_history, T, filename="velocity.mp4", interval=50, skip=3):
    fig, ax = plt.subplots()
    N = velocity_history[0][0].shape[0]
    x = np.linspace(0, 2*np.pi, N)
    y = np.linspace(0, 2*np.pi, N)
    X, Y = np.meshgrid(x, y, indexing='ij')

    u0, v0 = velocity_history[0]
    Q = ax.quiver(X[::skip, ::skip], Y[::skip, ::skip],
                  u0[::skip, ::skip], v0[::skip, ::skip], scale=50)
    ax.set_xlim(0, 2*np.pi)
    ax.set_ylim(0, 2*np.pi)

    def update(frame):
        u, v = velocity_history[frame]
        Q.set_UVC(u[::skip, ::skip], v[::skip, ::skip])
        ax.set_title("Velocity field at t = {:.3f}".format(T[frame]))
        return Q,

    ani = animation.FuncAnimation(fig, update, frames=len(velocity_history),
                                  interval=interval, blit=True)

    ani.save(filename, writer='ffmpeg', dpi=150)
    plt.close()


N = 128
dx = 2*np.pi / N
dt = 0.05
nu = 0
steps = 60000

zeta_history = []
psi_history = []
velocity_history = []
T = []
snapshot_interval = 100

zeta = initialize_vorticity(N, kind='vortex')
zeta_history.append(zeta)
T.append(0.0)

for n in range(steps):
    zeta = rk4_step(zeta, dt, dx, nu)
    if n % snapshot_interval == 0:
        psi = solve_poisson(zeta, dx)
        u, v = compute_velocity(psi, dx)
        velocity_history.append((1e1*u, 1e1*v))
        T.append(n * dt)
        zeta_history.append(zeta.copy())
        psi_history.append(psi.copy())
        print(f"Step {n}")

create_animation(zeta_history, T, filename="vorticity.mp4", interval=50)
create_animation(psi_history, T, filename="stream_func.mp4", interval=50)
create_velocity_animation(velocity_history, T, filename="velocity.mp4", interval=50)
