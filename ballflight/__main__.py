# encoding: utf-8 or plist.
"""
ballflight-- simulate ball flight with varying initial conditions for testing purposes
"""

import sys
import os
import matplotlib
import numpy as np

matplotlib.use("TkAgg")

from argparse import ArgumentParser
from argparse import RawDescriptionHelpFormatter

from options import PitchOptions
from results import PitchResults

from constants import *

__all__ = []
__version__ = 1.0
__date__ = '2018-10-19'
__updated__ = '2018-10-19'


def direction_to_axis(spin_direction, spin_efficiency):
    if spin_efficiency == 0:
        return [0, 1, 0]

    hours, minutes = map(float, spin_direction.split(":"))
    percent = (hours * 60 + minutes) / (60 * 12)
    angle = (percent * np.pi * 2) + (np.pi / 2)

    x = np.sin(angle)
    z = -np.cos(angle)

    term = np.square(1 - spin_efficiency)
    numerator = term * (np.square(x) + np.square(z))
    y = np.sqrt(numerator / (1 - term))

    return np.array([x, y, z]) / np.linalg.norm([x, y, z])


def main(args):
    play_type = "BASEBALL" if not args.fastPitch else "FAST_PITCH_SOFTBALL"
    spin_efficiency = args.breakSpin / args.spinRate
    spin_axis = direction_to_axis(args.spinDirection, spin_efficiency)
    pitch_options = PitchOptions(
        play_type=play_type,
        release_speed=args.releaseSpeed,
        release_spin_rate=args.spinRate,
        release_spin_axis=spin_axis,
        release_launch_angle=args.launchAngle,
        release_heading_angle=args.headingAngle,
        pitch=True
    )

    pitch_result = PitchResults.throw_ball(pitch_options)

    if args.plot:
        PitchResults.plot_results(pitch_result, pitch_options)

    res = np.array(pitch_result.pitch_break) * M_TO_IN
    print "Horizontal Break (in): %f" % res[0]
    print "Vertical Break (in): %f" % res[1]


if __name__ == "__main__":

    argv = sys.argv

    program_name = os.path.basename(sys.argv[0])
    program_version = "v%s" % __version__
    program_build_date = str(__updated__)
    program_version_message = '%%(prog)s %s (%s)' % (program_version, program_build_date)
    program_shortdesc = __import__('__main__').__doc__.split("\n")[1]

    # Setup argument parser
    parser = ArgumentParser(description=program_shortdesc, formatter_class=RawDescriptionHelpFormatter)
    parser.add_argument('-V', '--version', action='version', version=program_version_message)

    parser.add_argument('-r', '--releaseSpeed',
                        help="The release speed of the pitch.",
                        metavar="SPEED", required=False, type=float)
    parser.add_argument('-s', '--spinRate',
                        help="The release spin rate of the pitch.",
                        metavar="SPINRATE", required=False, type=float)
    parser.add_argument('-d', '--spinDirection',
                        help="The spin direction of the pitch (12:00, for example).",
                        metavar="SPINDIR", required=False)
    parser.add_argument('-b', '--breakSpin',
                        help="The spin contributing to break for this pitch. Also known as effective spin.",
                        metavar="BREAKSPIN", required=False, type=float)
    parser.add_argument('-l', '--launchAngle',
                        help="The launch angle of the pitch (0 is horizontal).",
                        metavar="LAUNCHANGLE", required=False, type=float,
                        default=0.0)
    parser.add_argument('-a', '--headingAngle',
                        help="The heading angle of the pitch (0 is toward the catcher).",
                        metavar="HEADINGANGLE", required=False, type=float,
                        default=0.0)
    parser.add_argument('-f', '--fastPitch',
                        help="Indicates that the pitch was with a softball.",
                        required=False, default=False,
                        action="store_true")
    parser.add_argument('-p', '--plot',
                        help="Indicates the user would like the flight path plotted.",
                        required=False, default=False,
                        action="store_true")

    # Process arguments
    args = parser.parse_args()

    main(args)
