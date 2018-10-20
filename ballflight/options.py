# This file contains the Options objects for pitches and hits. These objects contain all the information about the
#   the conditions of the situation that allow for accurate calculations.

from constants import *
import numpy as np


class PitchOptions(object):
    def __init__(self, play_type, release_speed, release_spin_rate, release_spin_axis, release_launch_angle,
                 release_heading_angle, pitch=True):
        self.playType = play_type  # must be all caps "FAST_PITCH_SOFTBALL" or "SLOW_PITCH_SOFTBALL" or "BASEBALL"
        self.ballMass = None
        self.ballRadius = None
        self.plateDistance = None
        self.initialPositionZ = None
        self.releaseSpeed = release_speed  # in MPH
        self.releaseSpinRate = release_spin_rate  # in RPM
        self.releaseSpinAxis = release_spin_axis  # a 1x3 unit vector of the spin axis
        self.releaseLaunchAngle = release_launch_angle  # in degrees
        self.releaseHeadingAngle = release_heading_angle  # in degrees
        self.pitch = pitch  # Boolean that indicates whether the ball thrown was a pitch, the default being that it was

        self.ballMass, self.ballRadius, self.plateDistance, self.initialPositionZ = self._play_type_parameters()

    def _play_type_parameters(self):
        if self.playType == 'SLOW_PITCH_SOFTBALL' or \
                self.playType == 'FAST_PITCH_SOFTBALL':
            ball_mass = 0.191
            ball_radius = 0.07625
            plate_distance = 13.106  # 43 ft
            if self.pitch:
                initial_position_z = 0.9  # about 3 feet, chosen based on assuming release at slightly below hip level
            else:
                initial_position_z = 2.0  # about 6.5 feet
            # Future Iteration Development: Include an age parameter that adjusts the mass, radius, and plate distance.
        else:  # baseball is the default
            ball_mass = 5.0 * OZ_TO_KG
            ball_radius = 0.03683
            plate_distance = 18.44  # 60.5 ft
            if self.pitch:
                initial_position_z = 1.8  # about 6 feet
            else:
                initial_position_z = 2.4  # about 8 feet

        return ball_mass, ball_radius, plate_distance, initial_position_z


class SwingOptions(object):
    def __init__(self, play_type, bat_length, bat_weight, bat_speed, bat_angle, session_type, pitch_speed=None,
                 pitch_angle=None, pitch_spin_rate=None):
        self.playType = play_type  # must be all caps "FAST_PITCH_SOFTBALL" or "SLOW_PITCH_SOFTBALL" or "BASEBALL"
        self.batLength = bat_length * IN_TO_M  # input in inches
        self.batWeight = bat_weight * OZ_TO_KG  # input in ounces
        self.batSpeed = bat_speed * MPH_TO_MS  # input in MPH
        self.bat_angle = bat_angle * DEG_TO_RAD  # input in degrees
        self.rotation_matrix_bat_to_global = None  # a rotation matrix that translates global vectors into the bat frame
        self.MOI_at_ref = None
        self.ref_to_knob = None
        self.COM_to_knob = None
        # Future Iteration Development: The Bat class inside of ballflight has other optional parameters (MOI_at_ref,
        #                               ref_to_knob, COM_to_knob) that give more information about the bat's Moment of
        #                               Inertia Matrix. At this point in time, that information is not widely known, so
        #                               the inputs here are left as None. In the future, if that information became
        #                               available, they should be added as SwingOptions inputs and called into the Bat
        #                               class we create later, inside of SwingResults.
        self.batDiameter = None
        self.batCOR = None
        self.ballMass = None
        self.ballRadius = None
        self.sessionType = session_type  # must be all caps "PITCH" or "FRONT_TOSS" or "TEE"
        self.pitch_speed = pitch_speed  # input in MPH
        self.pitch_angle = pitch_angle  # input in degrees
        self.pitch_spin_rate = pitch_spin_rate  # input in RPM
        # The default for the incoming pitch parameters is None, in which case we assume fastball. If the user does
        #   want to specify these inputs, the model will be more correct. The user can specify whichever inputs he/she
        #   wants, for example if the incoming speed was provided but not spin rate, the user should input the speed,
        #   and the spin rate will be assumed.

        self.batDiameter, self.batCOR, self.ballMass, self.ballRadius = self._play_type_parameters()

        self.rotation_matrix_bat_to_global = self._get_rotation_matrix()

        # If one of the three incoming pitch parameters has not been specified, the function to assign them assumptions
        #   is called.
        self.pitch_speed, self.pitch_angle, self.pitch_spin_rate = self._incoming_pitch_parameters()

    def _play_type_parameters(self):
        """
        This function creates some default options about the bat and ball based upon the play type.
        :return: bat_diameter: [m]
        :return: bat_COR: the coefficient of restitution of the bat
        :return: ball_mass: [kg]
        :return: ball_radius: [m]
        """
        if self.playType == 'SLOW_PITCH_SOFTBALL' or \
                self.playType == 'FAST_PITCH_SOFTBALL':
            bat_diameter = 0.05715  # 2.25 inches
            bat_COR = 1.2
            ball_mass = 0.191
            ball_radius = 0.07625
            # Future Iteration Development: Include an age parameter that adjusts the mass and radius.
        else:  # default is baseball
            bat_diameter = 0.066675  # 2.625 inches
            bat_COR = 0.5
            ball_mass = 0.145
            ball_radius = 0.03683

        return bat_diameter, bat_COR, ball_mass, ball_radius

    def _get_rotation_matrix(self):
        """
        The point of this function is to create a rotation matrix for the based based upon its initial orientation.
        :return: rotation_matrix: a rotation matrix that translates vectors from the global frame to the bat frame
        """
        rotation_matrix_from_bat_angle = np.zeros([3, 3])

        angle = -(np.pi - self.bat_angle)
        # The input bag angle is given from the -x axis being 0 and the -y axis being -pi (see README for pictures and
        #   explanation), so in order to get the correct angle for the rotation matrix (given from the +x axis) we
        #   subtract it from pi.

        rotation_matrix_from_bat_angle[0][0] = np.cos(angle)
        rotation_matrix_from_bat_angle[0][2] = -1.0 * np.sin(angle)
        rotation_matrix_from_bat_angle[1][1] = 1
        rotation_matrix_from_bat_angle[2][0] = np.sin(angle)
        rotation_matrix_from_bat_angle[2][2] = np.cos(angle)

        return rotation_matrix_from_bat_angle

    def _incoming_pitch_parameters(self):
        """
        This function uses the session type and the play type to determine information about the default incoming pitch.
        :return: pitch_speed: a scalar [m/s]
        :return: pitch_angle: [deg]
        :return: pitch_spin_rate: [rad/s]
        """
        if self.sessionType == 'PITCH':
            # All live pitches also use a collegiate level approximation. in the future, this could be adjusted
            #       according to the user's input level of competition.
            if self.playType == 'FAST_PITCH_SOFTBALL':
                pitch_speed = -60.0 * MPH_TO_MS
                pitch_angle = 0.0
                pitch_spin_rate = 1000 / RAD_TO_RPM
            elif self.playType == 'SLOW_PITCH_SOFTBALL':
                pitch_speed = -15.0 * MPH_TO_MS
                pitch_angle = 30.0 * DEG_TO_RAD
                pitch_spin_rate = 500 / RAD_TO_RPM
            else:  # the default will be baseball
                pitch_speed = -85.0 * MPH_TO_MS
                pitch_angle = 6.0 * DEG_TO_RAD
                pitch_spin_rate = 2225 / RAD_TO_RPM

        elif self.sessionType == 'FRONT_TOSS':
            if self.playType == 'FAST_PITCH_SOFTBALL':
                pitch_speed = -25.0 * MPH_TO_MS
                pitch_angle = 0.0
                pitch_spin_rate = 650 / RAD_TO_RPM
            elif self.playType == 'SLOW_PITCH_SOFTBALL':
                pitch_speed = -15.0 * MPH_TO_MS
                pitch_angle = 30.0 * DEG_TO_RAD
                pitch_spin_rate = 500 / RAD_TO_RPM
            else:  # the default will be baseball
                pitch_speed = -33.0 * MPH_TO_MS
                pitch_angle = 6.0 * DEG_TO_RAD
                pitch_spin_rate = 800 / RAD_TO_RPM

        # Note that all the pitch speeds are negative because the ball is coming in towards the batter.
        #   This also means that the incoming pitch angles are not set as negative, because they will automatically be
        #       pointing towards the batter and down when the negative speed is taken into account.

        else:  # the default will be a tee or soft toss session
            pitch_speed = 0.0
            pitch_angle = 0.0
            pitch_spin_rate = 0.0

        # The following if loops check if the inputs have been specified. If they have, then they leave them as the
        #   user input, adjusting units only. If not, they set them as the model's assumptions.
        if self.pitch_speed is None:
            incoming_pitch_speed = pitch_speed
        else:
            incoming_pitch_speed = self.pitch_speed * MPH_TO_MS

        if self.pitch_angle is None:
            incoming_pitch_angle = pitch_angle
        else:
            incoming_pitch_angle = self.pitch_angle * DEG_TO_RAD

        if self.pitch_spin_rate is None:
            incoming_pitch_spin_rate = pitch_spin_rate
        else:
            incoming_pitch_spin_rate = self.pitch_spin_rate / RAD_TO_RPM

        return incoming_pitch_speed, incoming_pitch_angle, incoming_pitch_spin_rate
