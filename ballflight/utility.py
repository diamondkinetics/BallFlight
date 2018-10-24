# encoding: utf-8 or plist.
"""
ballflight-- simulate ball flight with varying initial conditions for testing purposes
"""

import matplotlib
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import numpy as np

matplotlib.use("TkAgg")

from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt

from options import PitchOptions
from results import PitchResults

from constants import *


def direction_to_axis(spin_direction, spin_efficiency):
    if spin_efficiency == 0:
        return [0, 1, 0]

    hours, minutes = map(float, spin_direction.split(":"))
    percent = (hours * 60 + minutes) / (60 * 12)
    angle = (percent * np.pi * 2) + (np.pi / 2)

    x = np.sin(angle)
    z = -np.cos(angle)

    # term = np.square(1 - spin_efficiency**2)
    # numerator = term * 1
    y = np.sqrt(np.square(1.0/spin_efficiency)-1.0)

    return np.array([x, y, z]) / np.linalg.norm([x, y, z])


def break_from_velo_and_spin(direction, velo, spin_rate):
    spin_axis = direction_to_axis(direction, 1)

    pitch_options = PitchOptions(
        play_type="BASEBALL",
        release_speed=velo,
        release_spin_rate=spin_rate,
        release_spin_axis=spin_axis,
        release_launch_angle=0,
        release_heading_angle=0,
        pitch=True
    )

    pitch_result = PitchResults.throw_ball(pitch_options)

    return np.array(pitch_result.pitch_break) * M_TO_IN


if __name__ == "__main__":
    direction = "6:30"
    velo_range = np.arange(70, 95, 1)
    spin_range = np.arange(1000, 3000, 100)
    Z = [break_from_velo_and_spin(direction, velo, spin)[1] for spin in np.nditer(spin_range) for velo in
         np.nditer(velo_range)]
    Z = np.array(Z).reshape((len(spin_range), len(velo_range)))
    velo_range, spin_range = np.meshgrid(velo_range, spin_range)

    fig = plt.figure()
    ax = fig.gca(projection='3d')
    ax.set_title("Vertical Break vs Velo and Spin Rate")
    ax.set_xlabel('Velocity (mph)')
    ax.set_ylabel('Spin Rate (RPM)')

    # Plot the surface.
    surf = ax.plot_surface(velo_range, spin_range, Z, cmap=cm.coolwarm,
                           linewidth=0, antialiased=False)

    # Customize the z axis.
    ax.zaxis.set_major_locator(LinearLocator(10))
    ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))

    # Add a color bar which maps values to colors.
    fig.colorbar(surf, shrink=0.5, aspect=5)

    plt.show()
