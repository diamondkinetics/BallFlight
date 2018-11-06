# This file takes in information about the ball and bat and gives all of the information about the flight of that ball,
#  whether it was hit or pitched.

from __future__ import division
from constants import *
import numpy as np
import numpy.linalg as la
from enum import Enum


class Ball(object):
    def __init__(self, mass, radius):
        self.mass = mass
        self.radius = radius
        self.area = np.pi * (radius ** 2)


class BallEnum(Enum):
    BASEBALL = Ball(mass=0.145, radius=0.03683)
    SOFTBALL = Ball(mass=0.191, radius=0.07625)
    SOFTBALL_11IN = Ball(mass=0.1701, radius=0.07425)


class Bat(object):
    def __init__(self, length, mass, diameter, COR, MOI_at_ref=None, ref_to_knob=None, COM_to_knob=None):
        self.length = length  # [m]
        self.mass = mass  # [kg]
        self.diameter = diameter  # [m]
        self.COR = COR

        # the parameters MOI_at_ref, ref_to_knob, and COM_to_knob -- if we do, we can perform a more accurate MOI
        # calculation. if not, we treat the bat as a uniform cylinder and assume the COM is 75% of the bat's length
        self.MOI_at_ref = MOI_at_ref
        self.ref_to_knob = ref_to_knob
        self.COM_to_knob = COM_to_knob

        self.MOI_at_COM, self.COM_to_knob = self._bat_MOI_info()

    def _bat_MOI_info(self):
        """
        this function provides information about the COM and the MOI matrix
        :return: a 3x3 MOI matrix [kgm2] and a distance [m] between the COM and the knob
        """
        if self.COM_to_knob is None:
            # if this information was not provided, we assume the COM to be at 75% of the length
            bat_COM_to_knob = 0.75 * self.length
        else:
            bat_COM_to_knob = self.COM_to_knob

        if self.MOI_at_ref is None or self.ref_to_knob is None:
            # without enough information, we must use a uniform cylinder calculation for moment of inertia
            bat_MOI_at_COM = np.zeros(9)
            bat_MOI_at_COM[0] = 0.5 * self.mass * (0.5 * self.diameter) ** 2  # X axis, the long axis
            bat_MOI_at_COM[4] = self.mass * (1.0 / 12.0 * (self.length ** 2) + (0.25 * (0.375 * self.diameter) ** 2))
            # YZ axis
            bat_MOI_at_COM[8] = self.mass * (1.0 / 12.0 * (self.length ** 2) + (0.25 * (0.375 * self.diameter) ** 2))
            # YZ axis

            # equations above for cylinder moment of inertia from Source 5 listed in the README

        else:  # if we have enough information (the MOI_at_red and ref_to_knob inputs are filled) we can use the
            #       parallel axis theorem for increased accuracy
            ref_to_COM = np.fabs(bat_COM_to_knob - self.ref_to_knob)

            bat_MOI_at_COM = np.zeros(9)
            bat_MOI_at_COM[0] = 0.5 * self.mass * (0.5 * self.diameter) ** 2  # X axis, the long axis
            bat_MOI_at_COM[4] = self.MOI_at_ref - (self.mass * ref_to_COM ** 2)  # YZ axis
            bat_MOI_at_COM[8] = self.MOI_at_ref - (self.mass * ref_to_COM ** 2)  # YZ axis

        return bat_MOI_at_COM, bat_COM_to_knob


class SwungBat(object):
    def __init__(self, bat, barrel_velocity, rotation_matrix_bat_to_global):
        self.bat = bat  # an object of Bat class that is being swung
        self.global_impact_velocity = barrel_velocity  # a 3x1 vector of the barrel velocity at impact [m/s]
        self.rotation_matrix_bat_to_global = rotation_matrix_bat_to_global  # a rotation matrix that transforms
        #                                                                       vectors from the bat frame to the
        #                                                                       global frame
        self.rotation_matrix_global_to_bat = np.transpose(rotation_matrix_bat_to_global)  # a rotation matrix that
        #                                                                       transforms vectors from the global frame
        #                                                                       to the bat frame
        self.global_impact_direction, self.bat_frame_impact_direction = self._swing_parameters()

    def _swing_parameters(self):
        """
        :return: global_dir_impact_bat: 3x1 unit vector of the direction of impact -- global frame
        :return: bat_frame_dir_impact_bat: 3x1 unit vector of the direction of impact -- bat frame
        :return: global_impact_velocity: velocity of bat upon impact (3x1 vector) [m/s] -- global frame
        """
        bat_frame_dir_impact_bat = np.matmul(self.rotation_matrix_global_to_bat, self.global_impact_velocity)
        bat_frame_dir_impact_bat[0] = 0  # done to ignore the tangential slip of the ball
        # Future Iteration Development: account for tangential movement in the swing -- this will affect velocity,
        #                               direction and spin. If we account for it correctly, it will increase accuracy
        #                               and specificity of results.
        bat_frame_dir_norm = np.linalg.norm(bat_frame_dir_impact_bat)
        bat_frame_dir_impact_bat = np.multiply(bat_frame_dir_impact_bat, 1 / bat_frame_dir_norm)

        global_dir_impact_bat = np.matmul(self.rotation_matrix_bat_to_global, bat_frame_dir_impact_bat)

        return global_dir_impact_bat, bat_frame_dir_impact_bat


class PitchedBall(Ball):
    def __init__(self, ball, velocity, spin_rate, spin_axis):
        super(PitchedBall, self).__init__(ball.mass, ball.radius)
        self.velocity = velocity  # 3x1 vector of the velocity [m/s] in the global frame
        self.spin_rate = spin_rate  # scalar number of the spin rate [rad/s]
        self.spin_axis = spin_axis  # 1x3 unit vector in the direction of the spin axis

    @classmethod
    def from_incoming_info(cls, ball, pitch_speed, pitch_angle, pitch_spin_rate):
        """
        The purpose of this function is to initialize a PitchedBall that is going to be hit.
        We need this function because we do not have the same information about an incoming pitch that we do for a
            tracked pitch.
        :param ball: an object of Ball class
        :param pitch_speed: incoming speed of pitch (scalar) [m/s]
        :param pitch_angle: incoming angle of pitch [rad]
        :param pitch_spin_rate: incoming spin rate [rad/s]
        :return: pitched_ball: an object of the PitchedBall class (corresponding to the pitch delivered to batter)
        """
        pitched_velocity = np.array([0, pitch_speed * np.cos(pitch_angle), pitch_speed * np.sin(pitch_angle)]).reshape(3, 1)

        pitched_ball = PitchedBall(ball, pitched_velocity, pitch_spin_rate, [0., 0., 0.])
        # We have spin axis set as [0, 0, 0] because that is not currently taken into account in our hit ball model,
        #       as we have no way of obtaining that information yet.
        # Future Iteration Development: add pitch_spin_axis as an input for from_incoming_info and work that
        #                               information into the impact model (so the incoming spin axis would
        #                               affect the outgoing spin).

        return pitched_ball


class BallFlightInitialConditions(PitchedBall):
    def __init__(self, ball, velocity, spin_rate, spin_axis, initial_position):
        super(BallFlightInitialConditions, self).__init__(ball, velocity, spin_rate, spin_axis)
        self.ball = ball  # had to include this line because without it, get_no_spin was throwing an error, not sure why
        self.initial_position = initial_position  # 1x3 vector of starting position [m]

    @classmethod
    def from_pitched_ball(cls, pitched_ball, initial_position):
        return BallFlightInitialConditions(pitched_ball, pitched_ball.velocity, pitched_ball.spin_rate,
                                           pitched_ball.spin_axis, initial_position)

    @classmethod
    def from_impact(cls, pitched_ball, swung_bat, offset, initial_position):
        """
        The purpose of this function is to initialize a BallFlightInitialConditions object for a hit ball.
        :param pitched_ball: an object of the PitchedBall class
        :param swung_bat: an object of the SwungBat class
        :param offset: the distance between the center of mass of the bat and ball on impact [m]
        :param initial_position: 3x1 array indicating the starting position of the ball (point of impact) [m]
        :return: ball_flight_initial_conditions: an object of the BallFLightInitialConditions class
        """
        velocity = cls._batted_ball_velocity(pitched_ball, swung_bat)
        spin_rate, spin_axis = cls._batted_ball_spin(pitched_ball, swung_bat, offset)

        ball_flight_initial_conditions = BallFlightInitialConditions(pitched_ball, velocity, spin_rate, spin_axis,
                                                                     initial_position)

        return ball_flight_initial_conditions

    # The following 3 functions are helper functions for the from_impact function,
    #           a classmethod that initializes BallFlightInitialConditions during a swing.

    @classmethod
    def _batted_ball_spin(cls, pitched_ball, swung_bat, COM_separation):
            """
            This is a helper function for make_impact.
            :param pitched_ball: ball of PitchedBall class
            :param swung_bat: bat of SwungBat class
            :param COM_separation: distance between center of mass of ball and center of mass of bat [m]
            :return: exit_spin_rate and exit_spin_axis: spin rate [rad/s] and spin axis (1x3 unit vec) resulting from the impact
            """
            # ======== Spin Rate Calculation ========
            first_term = ((0.4 - swung_bat.bat.COR) / 1.4) * pitched_ball.spin_rate
            second_coeff = ((1 + swung_bat.bat.COR) / 1.4) / pitched_ball.radius
            second_term = second_coeff * pitched_ball.velocity[1]
            #                                           seeking the normal component of the incoming velocity
            third_coeff = second_coeff * (COM_separation / pitched_ball.radius)
            third_term = third_coeff * pitched_ball.velocity[0]
            #                                           seeking the tangential component of incoming velo

            exit_spin_rate = first_term + second_term - third_term

            # equations for this were primarily obtained from Source 3 listed in the README

            # ======== Spin Axis Calculation ========
            # Note: we assume that the axis is perfectly lined up with the bat, as we don't have the ability to make it
            #       more accurate than this. This loop just serves to differentiate between top spin and back spin.
            if COM_separation >= 0:  # if the ball strikes the top half of the bat -> back spin
                exit_spin_axis = [-1, 0, 0]
            else:  # if ball strikes bottom half of the bat -> top spin
                exit_spin_axis = [1, 0, 0]
            # Future Iteration Development: Incorporate tangential slip in order to accurately determine the
            #                               exiting ball's spin axis.

            return exit_spin_rate, exit_spin_axis

    @classmethod
    def _batted_ball_velocity(cls, pitched_ball, swung_bat):
        """
        This is a helper function for make_impact that finds the exit velocity of the ball after impact (in the global
            frame).
        :param pitched_ball: object of PitchedBall class
        :param swung_bat: object of SwungBat class
        :return: exit_velocity: a 3x1 vector of the ball's velocity [m/s]
        """
        bat_frame_ball_speed_norm, bat_frame_bat_speed_norm = cls._global_to_bat_frame(pitched_ball, swung_bat)

        B = 0.8 * swung_bat.bat.length - swung_bat.bat.COM_to_knob
        mass_ratio = (pitched_ball.mass / swung_bat.bat.mass) + \
                     (pitched_ball.mass * B ** 2 / swung_bat.bat.MOI_at_COM[8])
        den = 1.0 + mass_ratio
        num = bat_frame_ball_speed_norm * (mass_ratio - swung_bat.bat.COR) + (1 + swung_bat.bat.COR) *\
            bat_frame_bat_speed_norm

        exit_speed = num / den

        exit_velocity = np.multiply(exit_speed, swung_bat.global_impact_direction)

        return exit_velocity

    @classmethod
    def _global_to_bat_frame(cls, pitched_ball, swung_bat):
        """
        This is a helper function for batted_ball_velocity. It takes the 3D global frame vectors and puts them into the
            bat frame in 2D so that we can obtain the batted ball exit speed.
        :param pitched_ball: object of PitchedBall class
        :param swung_bat: object of SwungBat class
        :return: bat_frame_ball_speed_norm: the speed (scalar) of the bat in the bat frame in direction of impact [m/s]
        :return: bat_frame_bat_speed_norm: the speed (scalar) of the ball in the bat frame in direction of impact [m/s]
        """
        ball_velocity = pitched_ball.velocity

        bat_frame_ball_velocity = np.matmul(swung_bat.rotation_matrix_global_to_bat, ball_velocity)
        bat_frame_bat_velocity = np.matmul(swung_bat.rotation_matrix_global_to_bat, swung_bat.global_impact_velocity)

        bat_frame_ball_speed_norm = np.dot(np.transpose(bat_frame_ball_velocity), swung_bat.bat_frame_impact_direction)
        bat_frame_bat_speed_norm = np.dot(np.transpose(bat_frame_bat_velocity), swung_bat.bat_frame_impact_direction)

        return bat_frame_ball_speed_norm, bat_frame_bat_speed_norm

    def process(self, interval):
        """
        The purpose of this function is to initialize a BallFlightResults object.
        :param interval: time interval of measurements (in seconds)
        :return: flight_results: BallFlightResults object
        """
        # establish initial conditions
        dt = interval
        mass = self.mass
        radius = self.radius
        area = self.area
        spin_rate = self.spin_rate
        spin_axis = self.spin_axis
        pos_init = self.initial_position

        v_init = self.velocity
        grav = np.array([[0.], [0.], [-GRAVITY]])
        a_drag_init = self._acceleration_air_drag(v_init, mass, area)
        a_mag_init = self._magnus_accel(v_init, spin_rate, spin_axis, mass, area, radius)

        # motion of ball -- vector notation
        a_init_nograv = np.add(a_drag_init, a_mag_init)
        a_init = np.add(a_init_nograv, grav)

        pos_prev = pos_init
        v_prev = v_init
        a_prev = a_init

        pos = [pos_init]  # initializing empty array for position and velocity vectors
        vel = [v_init]

        # model flight
        while pos_prev[2] >= 0:  # loop terminates when ball hits the ground
            square_coeff = 0.5 * np.square(dt)
            pos_square = np.multiply(a_prev, square_coeff)
            pos_linear = np.multiply(v_prev, dt)
            pos_squ_lin = np.add(pos_square, pos_linear)
            pos_cur = np.add(pos_squ_lin, pos_prev)
            v_cur = np.add(np.multiply(a_prev, dt), v_prev)

            # update acceleration for next iterations
            a_drag_cur = self._acceleration_air_drag(v_cur, mass, area)
            a_mag_cur = self._magnus_accel(v_cur, spin_rate, spin_axis, mass, area, radius)
            a_cur_nograv = np.add(a_drag_cur, a_mag_cur)
            a_cur = np.add(a_cur_nograv, grav)

            # update previous pos and velo for next iteration

            pos_prev = pos_cur
            v_prev = v_cur
            a_prev = a_cur

            # append array
            pos.append(pos_cur)
            vel.append(v_cur)

        ball_flight_results = BallFlightResults(self, pos, vel)

        return ball_flight_results

    def process_no_spin(self, interval):
        """
        The purpose of this function is to initialize a BallFlightResults class from the same pitch, just without
            any of the spin data.
        :param interval: time interval of measurements [s]
        :return: no_spin_flight_results: an object of BallFlightResults class
        """
        no_spin_ball_flight_initial_conditions = self.get_no_spin()

        no_spin_ball_flight_results = no_spin_ball_flight_initial_conditions.process(interval)

        return no_spin_ball_flight_results

    def get_no_spin(self):
        """
        :return: a BallFlightInitialConditions object that has the same properties as the input, except it has no spin
        """
        ball_flight_initial_conditions_no_spin = BallFlightInitialConditions(self.ball, self.velocity, 0, [0, 0, 0],
                                                                             self.initial_position)
        return ball_flight_initial_conditions_no_spin

    # The following 3 functions are all helper functions for the process function,
    #           a function that initializes BallFlightResults.

    @classmethod
    def _acceleration_air_drag(cls, velocity, mass, area):
        """
        This is a helper function for flight_model.
        :param velocity: ball's velocity (1x3 vector) [m/s]
        :param mass: ball's mass [kg]
        :param area: ball's cross sectional area [m^2]
        :return: air_drag_accel: acceleration due to air drag (1x3 vector)
        """
        drag_magnitude = -0.5 * AIR_DENSITY * area * la.norm(velocity) * DRAG_COEFF
        air_drag = np.multiply(drag_magnitude, velocity)
        air_drag_accel = np.divide(air_drag, mass)

        # Note: the division of mass here is arbitrary. we use acceleration in our calculations, so we need to convert
        #       from force to acceleration at some point, but the choice to do it here is arbitrary.

        return air_drag_accel

    @classmethod
    def _coefficient_of_lift(cls, velocity, eff_spin_rate, radius, method ='NATHON'):
        """
        This is a helper function for flight_model.
        :param velocity: ball's velocity (3x1 vector) [m/s]
        :param eff_spin_rate: scalar value of ball's spin rate [rad/s]
        :param radius: ball's radius [m]
        :method: 'NATHON' or 'SAWICKI': different way for CL calculation
        :return: coefficient_of_lift: a scalar value
        """
        spin_parameter_num = radius * eff_spin_rate
        spin_parameter_den = la.norm(velocity)
        if spin_parameter_den == 0:
            spin_parameter = 0
        else:
            spin_parameter = spin_parameter_num / spin_parameter_den

        if method == 'NATHON':
            coeff_of_lift = spin_parameter/(0.4 + 2.32*spin_parameter)
        elif method == 'SAWICKI':
            if spin_parameter <= 0.1:
                coeff_of_lift = 1.5 * spin_parameter
            else:
                coeff_of_lift = 0.09 + 0.6 * spin_parameter

        # equation obtained from Source 2 listed in the README

        return coeff_of_lift

    @classmethod
    def _magnus_accel(cls, velocity, spin_rate, spin_axis, mass, area, radius):
        """
        This is a helper function for flight_model.
        :param velocity: ball's velocity (3x1 vector) [m/s]
        :param spin_rate: scalar value of ball's spin rate [rpm]
        :param spin_axis: direction of spin (1x3 unit vector)
        :param mass: ball's mass [kg]
        :param area: ball's cross sectional area [m^2]
        :param radius: ball's radius [m]
        :returns: magnus_accel: a 3x1 vector of the magnus acceleration
        """
        if np.all(spin_axis == 0) or spin_rate == 0 or np.all(velocity == 0):
            # want to avoid dividing by zero in the case of no spin
            magn_accel = np.array([[0.], [0.], [0.]])
        else:
            # get transverse component of spin
            spin_axis_parallel = np.dot(spin_axis, velocity.reshape(3)) * velocity.reshape(3) / (np.linalg.norm(velocity)**2)
            spin_axis_transverse = spin_axis - spin_axis_parallel
            eff_spin_rate = np.linalg.norm(spin_axis_transverse)*spin_rate

            omega = np.multiply(spin_rate, spin_axis)
            # creating a vector with magnitude of spin rate in direction of axis
            coeff_of_lift = cls._coefficient_of_lift(velocity, eff_spin_rate, radius, method='NATHON')
            magnus_mag = coeff_of_lift * 0.5 * AIR_DENSITY * area * np.asscalar(np.matmul(np.transpose(velocity),
                                                                                          velocity))
            # magnitude of the magnus force in each direction
            cross = np.cross(omega, np.transpose(velocity))
            cross_mag = la.norm(cross)
            unit_vec = np.multiply(cross, 1 / cross_mag)  # creating a unit vector for the direction of the magnus force
            magnus_force = magnus_mag * unit_vec
            magn_accel = np.transpose(np.multiply(magnus_force, 1 / mass))

        # equations for this were obtained from a variety of sources, listed in the README

        return magn_accel


class BallFlightResults(BallFlightInitialConditions):
    def __init__(self, ball_flight_initial_conditions, flight_position, flight_velocity):
        super(BallFlightResults, self).__init__(ball_flight_initial_conditions, ball_flight_initial_conditions.velocity,
                                                ball_flight_initial_conditions.spin_rate,
                                                ball_flight_initial_conditions.spin_axis,
                                                ball_flight_initial_conditions.initial_position)

        self.flight_position = flight_position  # this is an array of all of the position vectors throughout flight
        self.flight_velocity = flight_velocity  # this is an array of all the velocity vectors throughout flight

    def get_flight_over_distance(self, distance):
        """
        The flight model calculates the path of the ball until it hits the ground. Some pitches hit before the
            plate, and if that is the case, the last entry in the position array is good for the "plate" position of
            the ball, used to calculate break. However, if it does not hit before the plate, we need the entry that
            corresponds with the ball crossing the plate in order to calculate break. This function serves to terminate
            the position and velocity arrays at the point where the ball crosses the strike zone.

        :param distance: distance at which you'd like to stop the position/vector arrays [m]
                            - often would be the distance from the front of the plate plate to the pitchers mound
        :return: ball_flight_results: a BallFlightResults object with position and velocity vectors of the
                                        appropriate size
        """
        final_position = self.flight_position[-1]
        if final_position[1] < distance:
            pos_finals = self.flight_position
            vel_finals = self.flight_velocity
        else:
            index = len(self.flight_position)  # set the index as the last index so that if something goes wrong with
            #                                       the above portion of the loop, we still get the full position output
            for i in range(len(self.flight_position)):
                if self.flight_position[i][1] > distance:
                    index = i
                    break
            pos_finals = self.flight_position[:index + 1]
            vel_finals = self.flight_velocity[:index + 1]

        ball_flight_results = BallFlightResults(self, pos_finals, vel_finals)

        return ball_flight_results

