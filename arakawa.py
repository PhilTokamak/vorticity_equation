import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

def laplacian(f, dx):
    return (
        np.roll(f, 1, axis=0) + np.roll(f, -1, axis=0) +
        np.roll(f, 1, axis=1) + np.roll(f, -1, axis=1) - 4*f
    ) / dx**2

def arakawa_jacobian(z, psi, dx):
    # Arakawa Jacobian J_A = (J1 + J2 + J3) / 3
    J1 = (
        (np.roll(z, -1, axis=0) - np.roll(z, 1, axis=0)) *
        (np.roll(psi, -1, axis=1) - np.roll(psi, 1, axis=1)) -
        (np.roll(z, -1, axis=1) - np.roll(z, 1, axis=1)) *
        (np.roll(psi, -1, axis=0) - np.roll(psi, 1, axis=0))
    )

    J2 = (
        np.roll(z, -1, axis=0) * (
            np.roll(psi, -1, (0, 1)) - np.roll(psi, 1, (0, 1))
        ) -
        np.roll(z, 1, axis=0) * (
            np.roll(psi, -1, (0, 1)) - np.roll(psi, 1, (0, 1))
        ) -
        np.roll(z, -1, axis=1) * (
            np.roll(psi, -1, (1, 0)) - np.roll(psi, 1, (1, 0))
        ) +
        np.roll(z, 1, axis=1) * (
            np.roll(psi, -1, (1, 0)) - np.roll(psi, 1, (1, 0))
        )
    )

    J3 = (
        (np.roll(z, (-1, -1), (0, 1)) - np.roll(z, (-1, 1), (0, 1)) -
         np.roll(z, (1, -1), (0, 1)) + np.roll(z, (1, 1), (0, 1))) *
        (np.roll(psi, (-1, -1), (0, 1)) - np.roll(psi, (-1, 1), (0, 1)) -
         np.roll(psi, (1, -1), (0, 1)) + np.roll(psi, (1, 1), (0, 1)))
    )

    return (J1 + J2 + J3) / (12 * dx**2)

def solve_poisson(zeta, dx, tol=1e-6, max_iter=10000):
    psi = np.zeros_like(zeta)
    for _ in range(max_iter):
        psi_new = 0.25 * (
            np.roll(psi, 1, axis=0) + np.roll(psi, -1, axis=0) +
            np.roll(psi, 1, axis=1) + np.roll(psi, -1, axis=1) +
            dx**2 * zeta
        )
        if np.max(np.abs(psi_new - psi)) < tol:
            break
        psi = psi_new
    return psi

def rhs(zeta, dx, nu):
    psi = solve_poisson(-zeta, dx)
    jac = arakawa_jacobian(zeta, psi, dx)
    diff = nu * laplacian(zeta, dx)
    return -jac + diff

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
        return np.exp(-r2 / 0.1)
    elif kind == 'random':
        return np.random.randn(N, N) * 0.1
    else:
        raise ValueError("Unknown kind")


def create_animation(zeta_history, dx, filename="vorticity.mp4", interval=50):
    fig, ax = plt.subplots()
    im = ax.imshow(zeta_history[0], origin='lower', cmap='seismic',
                   extent=[0, 2*np.pi, 0, 2*np.pi], animated=True)
    cbar = plt.colorbar(im, ax=ax)

    def update(frame):
        im.set_array(zeta_history[frame])
        ax.set_title(f"Step {frame}")
        return [im]

    ani = animation.FuncAnimation(fig, update, frames=len(zeta_history),
                                  interval=interval, blit=True)

    ani.save(filename, writer='ffmpeg', dpi=150)
    plt.close()


N = 128
dx = 2*np.pi / N
dt = 0.01
nu = 0
steps = 1000

zeta_history = []
snapshot_interval = 10

zeta = initialize_vorticity(N, kind='vortex')

for n in range(steps):
    zeta = rk4_step(zeta, dt, dx, nu)
    if n % snapshot_interval == 0:
        zeta_history.append(zeta.copy())
        print(f"Step {n}")

create_animation(zeta_history, dx, filename="vorticity.mp4", interval=50)