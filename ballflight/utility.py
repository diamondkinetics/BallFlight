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
        release_launch_angle=3,
        release_heading_angle=0,
        pitch=True
    )

    pitch_result = PitchResults.throw_ball(pitch_options)

    return np.array(pitch_result.pitch_break) * M_TO_IN


def write_to_file(direction, velo_range, spin_range, horiz_breaks, vertical_breaks):
    with open(direction.replace(":", "") + ".csv", "w") as f:
        f.write("Spin Direction @ " + direction + ", Spin Rate\n")
        spin_list = "," + ",,".join(map(str, spin_range[0]))
        horiz_vert = "Release Velocity," + ("(H),(V)," * len(velo_range[0]))
        f.write(spin_list + "\n")
        f.write(horiz_vert + " \n")
        for i in range(velo_range.shape[0]):
            f.write(str(velo_range[i][0]) + ",")
            for j in range(velo_range.shape[1]):
                f.write(str(horiz_breaks[i][j]) + "," + str(vertical_breaks[i][j]) + ",")
            f.write("\n")


def for_direction(direction, plot):
    velo_range = np.arange(60, 95, 1)
    spin_range = np.arange(1000, 3000, 100)
    all_breaks = [break_from_velo_and_spin(direction, velo, spin) for velo in np.nditer(velo_range) for spin in np.nditer(spin_range)]
    horiz_breaks, vertical_breaks = zip(*all_breaks)
    spin_range, velo_range = np.meshgrid(spin_range, velo_range)
    vertical_breaks = np.array(vertical_breaks).reshape(velo_range.shape)
    horiz_breaks = np.array(horiz_breaks).reshape(velo_range.shape)

    write_to_file(direction, velo_range, spin_range, horiz_breaks, vertical_breaks)

    if plot:
        fig = plt.figure()
        ax = fig.gca(projection='3d')
        ax.set_title("Vertical Break vs Velo and Spin Rate")
        ax.set_xlabel('Velocity (mph)')
        ax.set_ylabel('Spin Rate (RPM)')

        # Plot the surface.
        surf = ax.plot_surface(velo_range, spin_range, vertical_breaks, cmap=cm.coolwarm,
                               linewidth=0, antialiased=False)

        # Customize the z axis.
        ax.zaxis.set_major_locator(LinearLocator(10))
        ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))

        # Add a color bar which maps values to colors.
        fig.colorbar(surf, shrink=0.5, aspect=5)

        plt.show()


if __name__ == "__main__":
    for hour in range(1, 13):
        for minutes in ["00", "30"]:
            print str(hour)+":"+minutes
            for_direction(str(hour)+":"+minutes, False)
