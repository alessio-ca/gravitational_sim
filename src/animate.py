import numpy as np
import matplotlib.pyplot as plt


class Animate:
    def __init__(self, ax, system, n_steps):
        self.ax = ax
        self.system = system
        self.n_steps = n_steps

        self.history_x = np.zeros(shape=(self.n_steps, self.system.n))
        self.history_y = np.zeros(shape=(self.n_steps, self.system.n))
        self.history_z = np.zeros(shape=(self.n_steps, self.system.n))

        prop_cycle = plt.rcParams["axes.prop_cycle"]
        colors = prop_cycle.by_key()["color"]

        self.plots = [
            self.ax.plot(
                planet.p[0], planet.p[1], planet.p[2], marker="o", color=colors[i]
            )[0]
            for i, planet in enumerate(self.system.planets)
        ]
        self.lines = [
            self.ax.plot([], [], [], color=colors[i])[0] for i in range(self.system.n)
        ]

    def __call__(self, j):
        # Perform motion (1 step)
        self.system.motion(n_steps=10)

        # Update
        for i, planet in enumerate(self.system.planets):
            # Append to history
            self.history_x[j, i] = planet.p[0]
            self.history_y[j, i] = planet.p[1]
            self.history_z[j, i] = planet.p[2]
            # Update lines
            self.lines[i].set_data(
                self.history_x[: (j + 1), i], self.history_y[: (j + 1), i]
            )
            self.lines[i].set_3d_properties(self.history_z[: (j + 1), i])
            # Update scatter
            self.plots[i].set_data(planet.p[0], planet.p[1])
            self.plots[i].set_3d_properties(planet.p[2])

        return (
            self.lines,
            self.plots,
        )
