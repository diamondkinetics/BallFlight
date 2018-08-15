# This file provides some plotting functions that will help to visualize the results the model returns.
from matplotlib import pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D


def plot_3d_ball_path(ball_results, ball_labels, viewType=None, ball_no_spin_results=None, pitch=False, softball=False):
    """
    :param ball_results: an array of all of the balls to be plotted (all of class BallFlightResults)
    :param ball_labels: an array of the corresponding labels
    :param viewType: a string (either 'TOP', 'SIDE', 'CATCHER', or 'PITCHER')
    :param ball_no_spin_results: an array of corresponding balls with no spin
    :param pitch: a Boolean that indicates if the ball being plotted is a pitch
    :param softball: a Boolean indicating either baseball or softball which allows for the appropriate field or
                        strike zone to be plotted
    :return: a plot of the ball path
    """
    fig = plt.figure(1)
    ax = fig.gca(projection='3d')

    colors = ['b-', 'g-', 'r-', 'c-', 'm-', 'y-', 'k-']  # this is arbitrarily ordered, just to provide colors for the
    #                                                       following loop. if more than 7 pitches are desired, this
    #                                                       array could be extended.

    for i in range(len(ball_results)):
        path = ball_results[i].flight_position
        path_x = grab_data(path, 0)
        path_y = grab_data(path, 1)
        path_z = grab_data(path, 2)
        if pitch:
            ax.plot(path_y, path_x, path_z, colors[i], label=ball_labels[i])
        else:
            ax.plot(path_x, path_y, path_z, colors[i], label=ball_labels[i])

        if ball_no_spin_results:
            no_spin_path = ball_no_spin_results[i].flight_position
            no_spin_path_x = grab_data(no_spin_path, 0)
            no_spin_path_y = grab_data(no_spin_path, 1)
            no_spin_path_z = grab_data(no_spin_path, 2)
            if pitch:
                ball_break = [round((path_x[-1] - no_spin_path_x[-1]) * 39.3701, 2), round((path_z[-1] -
                                                                                            no_spin_path_z[-1]) *
                                                                                           39.3701, 2)]
            else:
                ball_break = [round((path_x[-1] - no_spin_path_x[-1]) * 39.3701, 2), round((path_y[-1] -
                                                                                            no_spin_path_y[-1]) *
                                                                                           39.3701, 2)]
            ax.plot(no_spin_path_y, no_spin_path_x, no_spin_path_z, colors[i] + '-', label=ball_labels[i] +
                    ' No Spin - Break ' + str(ball_break))
    ax.legend()

    if viewType == 'TOP':
        ax.view_init(azim=-90, elev=90)
        ax.text2D(0.05, 0.95, "Top View", transform=ax.transAxes)
    elif viewType == 'SIDE':
        ax.view_init(azim=-90, elev=0)
        ax.text2D(0.05, 0.95, "Side View", transform=ax.transAxes)
    elif viewType == 'CATCHER':
        ax.view_init(azim=0, elev=0)
        ax.text2D(0.05, 0.95, "Catcher's View", transform=ax.transAxes)
    elif viewType == 'PITCHER':
        ax.view_init(azim=-180, elev=0)
        ax.text2D(0.05, 0.95, "Pitcher's View", transform=ax.transAxes)

    if softball:
        if pitch:
            plot_3d_softball_strike_zone(fig)
        else:  # default is not a pitch
            plot_3d_softball_field(fig)
    else:  # default is baseball
        if pitch:
            plot_3d_baseball_strike_zone(fig)
        else:
            plot_3d_baseball_field(fig)

    plt.show()
    return fig


def plot_pitches_in_zone(ball_results, ball_labels, plate_distance, ball_no_spin_results=None, softball=False):
    """
    :param ball_results: an array of all of the balls to be plotted (all of class BallFlightResults)
    :param ball_labels: an array of the corresponding labels
    :param plate_distance: the distance between mound and plate [m]
    :param ball_no_spin_results: an array of corresponding balls with no spin
    :param softball: a Boolean that indicates whether the sport is softball or not, which helps determine the
                        strikezone location
    :return: a plot of the ball path
    """
    plt.figure(2)
    colors = ['b.', 'g.', 'r.', 'c.', 'm.', 'y.', 'k.']  # this is arbitrarily ordered, just to provide colors for the
    #                                                       following loop. if more than 7 pitches are desired, this
    #                                                       array could be extended. the no_spin_colors array just
    #                                                       matches the order in this array so as to provide some visual
    #                                                       clarity as to the association between pitches and their
    #                                                       respective no spins.
    no_spin_colors = ['xkcd:light blue', 'xkcd:light green', 'xkcd:light red', 'xkcd:light cyan',
                      'xkcd:pale magenta', 'xkcd:light yellow', 'xkcd:grey']

    flight_distance = plate_distance - 1.524 - 0.4318  # we subtract out an average extension (5ft, maybe could be
    #                                                       customizable by height?) and plate length because
    #                                                       plate-to-mound measurements are from the plate apex, and
    #                                                       we care about where the ball crosses the front plane of the
    #                                                       plate.

    for i in range(len(ball_results)):
        pitch = ball_results[i].get_flight_over_distance(flight_distance)
        final_pos = pitch.flight_position[-1]
        plt.plot(-final_pos[0], final_pos[2], colors[i], label=ball_labels[i], markersize=14)
        if ball_no_spin_results:
            no_spin_pitch = ball_no_spin_results[i].get_flight_over_distance(flight_distance)
            no_spin_final_pos = no_spin_pitch.flight_position[-1]

            ball_break = [-round((final_pos[0] - no_spin_final_pos[0]) * 39.3701, 2), round((final_pos[2] -
                                                                                            no_spin_final_pos[2]) *
                                                                                            39.3701, 2)]
            plt.plot(-no_spin_final_pos[0], no_spin_final_pos[2], no_spin_colors[i], label=ball_labels[i] +
                     ' No Spin - Break ' + str(ball_break), marker='.', markersize=14)
    if softball:
        # strike zone is as wide as the plate (17in), approximately 1.25 feet off the ground, and about 2 feet in
        #       vertical direction
        # horizontal lines of strike zone
        plt.plot([-0.2159, 0.2159], [0.3810, 0.3810], 'k-')
        plt.plot([-0.2159, 0.2159], [0.5873, 0.5873], 'k-')
        plt.plot([-0.2159, 0.2159], [0.79366, 0.79366], 'k-')
        plt.plot([-0.2159, 0.2159], [1.0, 1.0], 'k-')
        # vertical lines of strike zone
        plt.plot([-0.2159, -0.2159], [0.3810, 1.0], 'k-')
        plt.plot([-0.07197, -0.07197], [0.3810, 1.0], 'k-')
        plt.plot([0.07197, 0.07197], [0.3810, 1.0], 'k-')
        plt.plot([0.2159, 0.2159], [0.3810, 1.0], 'k-')
    else:  # baseball zone is the default
        # strike zone is as wide as the plate (17in), approximately 2 feet off the ground, and about 2 feet in vertical
        #       direction
        # horizontal lines of strike zone
        plt.plot([-0.2159, 0.2159], [0.63846, 0.63846], 'k-')
        plt.plot([-0.2159, 0.2159], [0.82242, 0.82242], 'k-')
        plt.plot([-0.2159, 0.2159], [1.00638, 1.00638], 'k-')
        plt.plot([-0.2159, 0.2159], [1.19034, 1.19034], 'k-')
        # vertical lines of strike zone
        plt.plot([-0.2159, -0.2159], [0.63846, 1.19034], 'k-')
        plt.plot([-0.07197, -0.07197], [0.63846, 1.19034], 'k-')
        plt.plot([0.07197, 0.07197], [0.63846, 1.19034], 'k-')
        plt.plot([0.2159, 0.2159], [0.63846, 1.19034], 'k-')

    plt.xlim([-1.5875, 1.5875])  # wide enough to include both batters boxes
    plt.ylim([-0.1, 1.8288])  # as tall as average mlb player (6ft)

    plt.title('Strike Zone View (Pitcher)')

    plt.legend()

    plt.show()


def plot_3d_baseball_strike_zone(fig):
    """
    This function adds a strike zone to a 3D plot (assuming metric axes).
    :param: fig: a figure that the strikezone will go on
    """
    ax = fig.gca(projection='3d')

    # strike zone is as wide as the plate (17in), approximately 2 feet off the ground, and 2 feet in vertical direction
    # strike zone is set at plate distance (60ft) minus 5 ft extension and 17 in plate

    ax.plot([16.332, 16.332], [-0.2159, 0.2159], [0.63846, 0.63846], 'k-')
    ax.plot([16.332, 16.332], [0.2159, 0.2159], [0.63846, 1.19034], 'k-')
    ax.plot([16.332, 16.332], [-0.2159, 0.2159], [1.19034, 1.19034], 'k-')
    ax.plot([16.332, 16.332], [-0.2159, -0.2159], [0.63846, 1.19034], 'k-')

    ax.plot([16.332, 16.332], [-0.2159, 0.2159], [0.82242, 0.82242], 'k-')
    ax.plot([16.332, 16.332], [-0.2159, 0.2159], [1.00638, 1.00638], 'k-')
    ax.plot([16.332, 16.332], [-0.07197, -0.07197], [0.63846, 1.19034], 'k-')
    ax.plot([16.332, 16.332], [0.07197, 0.07197], [0.63846, 1.19034], 'k-')

    # sets axes at a reasonable pitch view
    ax.set_xlim3d([-1.524, 17.1958])
    ax.set_ylim3d([-1.5875, 1.5875])
    ax.set_zlim3d([0, 1.8288])

    return


def plot_3d_baseball_field(fig):
    """
    This function adds a baseball field to a 3D plot (assuming metric axes).
    :param: fig: a figure that the field will go on
    """
    ax = fig.gca(projection='3d')
    # the following lines are creating the outfield edge
    outfield_val = -70.046
    outfield_x = [outfield_val]
    while outfield_val < 69.954:  # not sure why, but adding one gets it slightly below value so i account for that by
        #                               shortening the loop
        outfield_val += 1
        outfield_x.append(outfield_val)

    outfield_coeffs = [-1.05726306e-02, 0., 1.21920000e+02]
    # linear fit for the outfield curve using the points (-70.046, 70.046), (0, 121.92), (70.046, 70.046) which
    #   corresponds with foul poles at 325 ft and center at 400 ft
    outfield_y = []
    for i in outfield_x:
        y = outfield_coeffs[0] * (i ** 2) + outfield_coeffs[1] * i + outfield_coeffs[2]
        outfield_y.append(y)
    ax.plot(outfield_x, outfield_y, np.zeros(len(outfield_x)), 'g-')

    # the following lines are creating the infield edge
    infield_val = -27.502
    infield_x = [infield_val]
    while infield_val < 27.498:  # same issue as in the above loop
        infield_val += 1
        infield_x.append(infield_val)

    infield_coeffs = [-2.650431e-02, 0., 4.754880e+01]
    # linear fit for the infield curve using the points (-27.502, 27.502), (0, 47.5488), (27.502, 27.502) which
    #   corresponds with an arc centered at the rubber with a radius of 95 feet
    infield_y = []
    for i in infield_x:
        y = infield_coeffs[0] * (i ** 2) + infield_coeffs[1] * i + infield_coeffs[2]
        infield_y.append(y)
    ax.plot(infield_x, infield_y, np.zeros(len(infield_x)), 'xkcd:brown')

    # the following lines are plotting the foul lines
    ax.plot([0, -70.046], [0, 70.046], [0, 0], 'k-')
    ax.plot([0, 70.046], [0, 70.046], [0, 0], 'k-')
    ax.set_zlim3d([0, 54.864])
    # z axis limited to 180 ft because this is one of the highest HRs ever recorded (Hanley Ramirez, 2015)


def plot_3d_softball_strike_zone(fig):
    """
    This function adds a strike zone to a 3D plot (assuming metric axes).
    :param: fig: a figure that the strikezone will go on
    """
    ax = fig.gca(projection='3d')

    # strike zone is as wide as the plate (17in), approximately 1.25 feet off the ground, and about 2 feet in vertical
    #       direction
    # strike zone is set at plate distance (43ft) minus 5 ft extension and 17 in plate
    ax.plot([11.151, 11.151], [-0.2159, 0.2159], [0.3810, 0.3810], 'k-')
    ax.plot([11.151, 11.151], [0.2159, 0.2159], [0.3810, 1.0], 'k-')
    ax.plot([11.151, 11.151], [-0.2159, 0.2159], [1.0, 1.0], 'k-')
    ax.plot([11.151, 11.151], [-0.2159, -0.2159], [0.3810, 1.0], 'k-')

    ax.plot([11.151, 11.151], [-0.2159, 0.2159], [0.5873, 0.5873], 'k-')
    ax.plot([11.151, 11.151], [-0.2159, 0.2159], [0.79366, 0.79366], 'k-')
    ax.plot([11.151, 11.151], [-0.07197, -0.07197], [0.3820, 1.0], 'k-')
    ax.plot([11.151, 11.151], [0.07197, 0.07197], [0.3820, 1.0], 'k-')

    # sets axes at a reasonable pitch view
    ax.set_xlim3d([-1.524, 17.1958])
    ax.set_ylim3d([-1.5875, 1.5875])
    ax.set_zlim3d([0, 1.8288])

    return


def plot_3d_softball_field(fig):
    """
    This function adds a softball field to a 3D plot (assuming metric axes).
    :param: fig: a figure that the field will go on
    """
    ax = fig.gca(projection='3d')
    # the following lines are creating the outfield edge
    outfield_val = -43.105
    outfield_x = [outfield_val]
    while outfield_val < 42.0:  # not sure why, but adding one gets it slightly below value so i account for that by
        #                               shortening the loop
        outfield_val += 1
        outfield_x.append(outfield_val)

    outfield_coeffs = [-1.28904581e-02, 0., 6.70560000e+01]
    # linear fit for the outfield curve using the points (-43.105, 43.105), (0, 67.056), (43.105, 43.105) which
    #   corresponds with foul poles at 200 ft and center at 220 ft
    outfield_y = []
    for i in outfield_x:
        y = outfield_coeffs[0] * (i ** 2) + outfield_coeffs[1] * i + outfield_coeffs[2]
        outfield_y.append(y)
    ax.plot(outfield_x, outfield_y, np.zeros(len(outfield_x)), 'g-')

    # the following lines are creating the infield edge
    infield_val = -17.701
    infield_x = [infield_val]
    while infield_val < 17.0:  # same issue as in the above loop
        infield_val += 1
        infield_x.append(infield_val)

    infield_coeffs = [-0.04370345, 0., 31.3944]
    # linear fit for the outfield curve using the points (-17.701, 17.701), (0, 38.140), (17.701, 17.701) which
    #   corresponds with an arc centered at the rubber with a radius of 60 feet
    infield_y = []
    for i in infield_x:
        y = infield_coeffs[0] * (i ** 2) + infield_coeffs[1] * i + infield_coeffs[2]
        infield_y.append(y)
    ax.plot(infield_x, infield_y, np.zeros(len(infield_x)), 'xkcd:brown')

    # the following lines are plotting the foul lines
    ax.plot([0, -43.105], [0, 43.105], [0, 0], 'k-')
    ax.plot([0, 43.105], [0, 43.105], [0, 0], 'k-')
    ax.set_zlim3d([0, 54.864])
    # z axis limited to 180 ft because this is one of the highest HRs ever recorded (Hanley Ramirez, 2015)


def grab_data(array, index):
    """
    This is a helper function for parsing position output.
    :param array: any array that's items are lists
    :param index: desired index of the sublists
    :return: subset: a list of all of the values at the given index in each sublist
    """
    subset = []
    for i in array:
        subset.append(float(i[index]))
    return subset
