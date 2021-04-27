import numpy as np
import matplotlib.pyplot as plt
from astropy.time import Time
from matplotlib.offsetbox import OffsetImage, AnnotationBbox
from mpl_toolkits.mplot3d import proj3d
from src.utils import load_parameters

PARAMS = load_parameters()


class Animate:
    def __init__(self, ax, system, n_steps, names, subnames=None):
        self.ax = ax
        self.system = system
        self.n_steps = n_steps
        self.names = names
        if subnames:
            self.subnames = subnames
        else:
            self.subnames = self.names
        self.num_planets = len(self.subnames)
        self.planets_to_plot = [
            planet
            for name, planet in zip(self.names, self.system.planets)
            if name in self.subnames
        ]

        self.history_x = np.zeros(shape=(self.n_steps, self.num_planets))
        self.history_y = np.zeros(shape=(self.n_steps, self.num_planets))
        self.history_z = np.zeros(shape=(self.n_steps, self.num_planets))

        self.time = Time.now().jd

        prop_cycle = plt.rcParams["axes.prop_cycle"]
        colors = prop_cycle.by_key()["color"]

        self.lines = [
            self.ax.plot(
                [], [], [], color=colors[i % len(colors)], label=self.subnames[i]
            )[0]
            for i in range(self.num_planets)
        ]

        self.timestamp = self.ax.text(
            x=ax.get_xlim()[0],
            y=1.5 * ax.get_ylim()[1],
            z=1.5 * ax.get_zlim()[1],
            s="Date: ",
            color="k",
            # transform=self.ax.transAxes,
            fontsize="x-large",
        )

    def __call__(self, j):
        # Perform motion (10 step, corresponds to 1 day)
        self.system.motion(10)
        self.time += 1
        # Update
        for i, planet in enumerate(self.planets_to_plot):

            # Append to history
            self.history_x[j, i] = planet.p[0]
            self.history_y[j, i] = planet.p[1]
            self.history_z[j, i] = planet.p[2]
            # Update lines
            self.lines[i].set_data(
                self.history_x[: (j + 1), i], self.history_y[: (j + 1), i]
            )
            self.lines[i].set_3d_properties(self.history_z[: (j + 1), i])

            # Update date
            self.timestamp.set_text(
                "Date: "
                + Time(
                    Time(self.time, format="jd"), format="iso", out_subfmt="date"
                ).iso
            )

        return (self.lines, self.timestamp)


class Animate_Point(Animate):
    def __init__(self, ax, system, n_steps, names, subnames=None):
        super().__init__(
            ax=ax, system=system, n_steps=n_steps, names=names, subnames=subnames
        )
        prop_cycle = plt.rcParams["axes.prop_cycle"]
        colors = prop_cycle.by_key()["color"]
        self.plots = [
            self.ax.plot(
                planet.p[0],
                planet.p[1],
                planet.p[2],
                marker="o",
                color=colors[i % len(colors)],
            )[0]
            for i, planet in enumerate(self.planets_to_plot)
        ]

    def __call__(self, j):
        self.lines, self.timestamp = super().__call__(j)
        # Update scatter
        for i, planet in enumerate(self.planets_to_plot):
            self.plots[i].set_data(planet.p[0], planet.p[1])
            self.plots[i].set_3d_properties(planet.p[2])

        return (self.lines, self.plots, self.timestamp)


class Animate_Icon(Animate):
    def __init__(self, ax, system, n_steps, names, subnames=None):
        super().__init__(
            ax=ax, system=system, n_steps=n_steps, names=names, subnames=subnames
        )

        self.annotations = []

        zooms = {
            name: value
            for name, value in PARAMS["zoom"].items()
            if name in self.subnames
        }
        ids = {
            name: value for name, value in PARAMS["id"].items() if name in self.subnames
        }
        paths = [
            f'{PARAMS["icon_paths"]}{num:03d}_{name}.png'
            for num, name in zip(ids.values(), self.subnames)
        ]
        for planet, path, zoom in zip(self.planets_to_plot, paths, zooms.values()):
            x2, y2, _ = proj3d.proj_transform(
                planet.p[0], planet.p[1], planet.p[2], self.ax.get_proj()
            )
            ab = AnnotationBbox(
                OffsetImage(plt.imread(path), zoom=zoom), (x2, y2), frameon=False
            )
            self.ax.add_artist(ab)
            self.annotations.append(ab)

    def __call__(self, j):
        self.lines, self.timestamp = super().__call__(j)
        # Update annotations
        for i, planet in enumerate(self.planets_to_plot):
            # Update annotations
            x2, y2, _ = proj3d.proj_transform(
                planet.p[0], planet.p[1], planet.p[2], self.ax.get_proj()
            )
            self.annotations[i].xybox = (x2, y2)
            # Update date
            self.timestamp.set_text(
                "Date: "
                + Time(
                    Time(self.time, format="jd"), format="iso", out_subfmt="date"
                ).iso
            )

        return (self.lines, self.annotations, self.timestamp)


class AnimateList:
    def __init__(self, ax, system, delta_step, names, subnames=None):
        self.ax = ax
        self.system = system
        self.dt = system.dt
        self.delta_step = delta_step
        self.names = names
        if subnames:
            self.subnames = subnames
        else:
            self.subnames = self.names
        self.num_planets = len(self.subnames)
        self.index_to_plots = [
            i for i, name in enumerate(self.names) if name in self.subnames
        ]

        self.time = Time.now().jd

        prop_cycle = plt.rcParams["axes.prop_cycle"]
        colors = prop_cycle.by_key()["color"]

        self.lines = [
            self.ax.plot(
                [], [], [], color=colors[i % len(colors)], label=self.subnames[i]
            )[0]
            for i in range(self.num_planets)
        ]

        self.timestamp = self.ax.text(
            x=ax.get_xlim()[0],
            y=1.5 * ax.get_ylim()[1],
            z=1.5 * ax.get_zlim()[1],
            s="Date: ",
            color="k",
            # transform=self.ax.transAxes,
            fontsize="x-large",
        )

    def __call__(self, j, positions):
        self.time += self.delta_step * self.dt
        # Update
        for i, planet in enumerate(self.index_to_plots):
            # Update lines
            self.lines[i].set_data(
                positions[: (j + 1) * self.delta_step : self.delta_step, planet, 0],
                positions[: (j + 1) * self.delta_step : self.delta_step, planet, 1],
            )
            self.lines[i].set_3d_properties(
                positions[: (j + 1) * self.delta_step : self.delta_step, planet, 2]
            )

            # Update date
            self.timestamp.set_text(
                "Date: "
                + Time(
                    Time(self.time, format="jd"), format="iso", out_subfmt="date"
                ).iso
            )

        return (self.lines, self.timestamp)


class AnimateList_Point(AnimateList):
    def __init__(self, ax, system, delta_step, names, subnames=None):
        super().__init__(
            ax=ax, system=system, delta_step=delta_step, names=names, subnames=subnames
        )
        prop_cycle = plt.rcParams["axes.prop_cycle"]
        colors = prop_cycle.by_key()["color"]
        self.plots = [
            self.ax.plot(
                system.positions[planet, 0],
                system.positions[planet, 1],
                system.positions[planet, 2],
                marker="o",
                color=colors[i % len(colors)],
            )[0]
            for i, planet in enumerate(self.index_to_plots)
        ]

    def __call__(self, j, positions):
        self.lines, self.timestamp = super().__call__(j, positions)
        # Update scatter
        for i, planet in enumerate(self.index_to_plots):
            x_pos = positions[j * self.delta_step, planet, 0]
            y_pos = positions[j * self.delta_step, planet, 1]
            z_pos = positions[j * self.delta_step, planet, 2]
            self.plots[i].set_data(x_pos, y_pos)
            self.plots[i].set_3d_properties(z_pos)

        return (self.lines, self.plots, self.timestamp)
