from __future__ import annotations
import numpy as np
from numba.experimental import jitclass
from numba.types import List
from numba import float64  # import the types

# Define astronomical units
AU = 149597870700
T = 86400.0
M = 1.98847 * 10 ** 30

G = (6.67430 / 10 ** 11) * M * T ** 2 / AU ** 3


@jitclass
class Planet:
    m: float
    p: float64[:]
    v: float64[:]
    f: float64[:]

    def __init__(self, m: float):
        self.m = m
        self.p = np.zeros(shape=(3,), dtype=np.float64)
        self.v = np.zeros(shape=(3,), dtype=np.float64)
        self.f = np.zeros(shape=(3,), dtype=np.float64)

    def get_force(self, planet: Planet):
        delta = planet.p - self.p
        dot_prod = np.sum(delta * delta)
        energy = G * self.m * planet.m / np.sqrt(dot_prod)
        force = energy * delta / dot_prod

        return force, energy

    def update_pos(self, dt: float):
        return self.p + self.v * dt + 0.5 * (self.f / self.m) * (dt) * (dt)

    def update_vel(self, dt: float):
        return self.v + 0.5 * (self.f / self.m) * dt


@jitclass
class System:
    n: int
    dt: float
    masses: float64[:]
    positions: float64[:, :]
    velocities: float64[:, :]
    U: float
    K: float
    min_idx: int
    planets: List(Planet.class_type.instance_type)

    def __init__(
        self,
        n: int,
        dt: float,
        m_in: np.ndarray,
        p_in: np.ndarray,
        v_in: np.ndarray,
        sun_static: bool = False,
    ):
        # Perform initialization checks
        assert n == m_in.shape[0]
        assert n == p_in.shape[0]
        assert n == v_in.shape[0]
        self.n = n
        self.dt = dt
        self.masses = m_in.astype(np.float64)
        self.positions = p_in.astype(np.float64)
        self.velocities = v_in.astype(np.float64)

        self.U = 0
        self.K = 0

        if sun_static:
            self.min_idx = 1
        else:
            self.min_idx = 0

        # Perform initialization
        self.reset()

    def reset(self):
        # This is a dummy initialization for numba to infer typing of list
        self.planets = [Planet(self.masses[0])]
        self.planets.pop()

        # Declare planets
        for i in range(self.n):
            body = Planet(self.masses[i])
            body.p = self.positions[i]
            body.v = self.velocities[i]
            self.planets.append(body)

        # Initialize forces and energy
        self.apply_forces()
        # Calculate kinetic energy
        self.calculate_kin()

        return self

    def apply_forces(self):
        # Zero the forces and energies
        for i in range(self.n):
            self.planets[i].f = np.zeros(shape=(3,), dtype=np.float64)

        self.U = 0

        # Calculate forces
        for i in range(self.n - 1):
            for j in range(i + 1, self.n):
                force, energy = self.planets[i].get_force(self.planets[j])
                self.planets[i].f += force
                self.planets[j].f -= force

                # Gravitational energy contributions are negative
                self.U -= energy

        return self

    def calculate_kin(self):
        self.K = 0
        for i in range(self.n):
            self.K += (
                0.5 * self.planets[i].m * np.sum(self.planets[i].v * self.planets[i].v)
            )

    def integration(self):
        # Update positions and velocities
        for i in range(self.min_idx, self.n):
            self.planets[i].p = self.planets[i].update_pos(self.dt)
            self.planets[i].v = self.planets[i].update_vel(self.dt)

        # Update acceleration
        self.apply_forces()

        # Update velocities
        for i in range(self.min_idx, self.n):
            self.planets[i].v = self.planets[i].update_vel(self.dt)

        # Calculate kinetic
        self.calculate_kin()

        return self

    def motion(self, n_steps: int):
        for _ in range(n_steps):
            self.integration()
        return self

    def motion_dump(self, n_steps: int):
        positions = np.zeros(shape=(n_steps, self.n, 3), dtype=np.float64)
        for i in range(n_steps):
            self.integration()
            for j, planet in enumerate(self.planets):
                positions[i, j, :] = planet.p

        return positions
