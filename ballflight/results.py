# This file uses the ballflight model to pull actual results about the pitch or hit. This is the file that users should
#   call from when they seek their metric output.

from ballflight import *
from options import *
import plotting as bfp


class PitchResults(object):
    def __init__(self):
        self.position_of_pitch = []
        self.position_of_nospin_pitch = []
        self.pitch_break = [0.0, 0.0]
        # First index is horizontal, second is vertical. Break is given from the perspective of the thrower. For
        #   example, if the ball is pitched, positive break is towards the righty batter's box. If the ball is thrown,
        #   positive break is towards right field.

    @classmethod
    def throw_ball(cls, pitch_options):
        """
        This function takes in the basics of a thrown ball and gives back the final resulting model.
        :param: pitch_options: an object of the PitchOptions class
        :return: pitch_result: an object of the PitchResult class, containing information about the path of the ball
                                and the final break.
        """
        pitch_result = PitchResults()

        velocity_z_component = (pitch_options.releaseSpeed * MPH_TO_MS) * np.sin(pitch_options.releaseLaunchAngle
                                                                                 * DEG_TO_RAD)
        velocity_xy_component = (pitch_options.releaseSpeed * MPH_TO_MS) * np.cos(pitch_options.releaseLaunchAngle
                                                                                  * DEG_TO_RAD)
        velocity_x_component = velocity_xy_component * np.sin(pitch_options.releaseHeadingAngle * DEG_TO_RAD)
        velocity_y_component = velocity_xy_component * np.cos(pitch_options.releaseHeadingAngle * DEG_TO_RAD)

        velocity = np.array([[velocity_x_component], [velocity_y_component], [velocity_z_component]])

        ball = Ball(pitch_options.ballMass, pitch_options.ballRadius)

        ball_pitched = PitchedBall(ball, velocity, pitch_options.releaseSpinRate / RAD_TO_RPM,
                                   pitch_options.releaseSpinAxis)

        init_position = np.array([[0.], [0.], [pitch_options.initialPositionZ]])

        ball_flight_initial_conditions = BallFlightInitialConditions.from_pitched_ball(ball_pitched, init_position)

        spin_flight_results = ball_flight_initial_conditions.process(0.001)
        # Note: a time step interval of 0.001 seconds was chosen here arbitrarily as a time that provided a sufficient
        #       yet not excessive number of data points.

        no_spin_flight_results = ball_flight_initial_conditions.process_no_spin(0.001)

        # *********** Break ***********
        if pitch_options.pitch:  # if the throw has been specified as a pitch, we terminate the flight at the plate
            #                       distance and calculate break in the vertical and horizontal directions in terms of
            #                       the strike zone ([2, 3] would mean two inches right and 3 inches up.
            # flight_distance = pitch_options.plateDistance - 1.524 - 0.4318  # we subtract out an average extension
            flight_distance = 14
            #                                                                   (5ft) and plate length because
            #                                                                   plate-to-mound measurements are from
            #                                                                   the plate apex, and we care about where
            #                                                                   the ball crosses the front plane of the
            #                                                                   plate.

            spin_flight_results = spin_flight_results.get_flight_over_distance(flight_distance)
            no_spin_flight_results = no_spin_flight_results.get_flight_over_distance(flight_distance)
            final_position_spun = spin_flight_results.flight_position[-1]
            final_position_unspun = no_spin_flight_results.flight_position[-1]

            x_direction_break = final_position_spun[0]-final_position_unspun[0]
            z_direction_break = final_position_spun[2]-final_position_unspun[2]
            pitch_break = [-x_direction_break[0], z_direction_break[0]]

        else:  # if the throw is not a pitch, we allow the flight to continue until it hits the ground, at which point
            #        we report break in the x and y, where x is deviation from the center line (positive x would be
            #        towards right), and y is distance straight out from the plate.
            final_position_spun = spin_flight_results.flight_position[-1]
            final_position_unspun = no_spin_flight_results.flight_position[-1]
            x_direction_break = final_position_spun[0] - final_position_unspun[0]
            y_direction_break = final_position_spun[1] - final_position_unspun[1]
            pitch_break = [-x_direction_break[0], y_direction_break[0]]

        pitch_result.position_of_pitch = spin_flight_results
        pitch_result.position_of_nospin_pitch = no_spin_flight_results
        pitch_result.pitch_break = pitch_break

        return pitch_result


    @classmethod
    def plot_results(cls, pitch_result, pitch_options):
        # Plotting Functions
        if pitch_options.playType == 'SLOW_PITCH_SOFTBALL' or \
                        pitch_options.playType == 'FAST_PITCH_SOFTBALL':
            softball = True
        else:
            softball = False

        bfp.plot_3d_ball_path([pitch_result.position_of_pitch], [''], None, [pitch_result.position_of_nospin_pitch],
                              pitch_options.pitch, softball)

        if pitch_options.pitch:
            bfp.plot_pitches_in_zone([pitch_result.position_of_pitch], [''], pitch_options.plateDistance, [pitch_result.position_of_nospin_pitch],
                                     softball)


class SwingResults(object):
    def __init__(self):
        self.flight_position = []
        self.max_distance = 0.0
        self.ideal_exit_speed = 0.0
        self.ideal_launch_angle = 0.0
        self.ideal_heading_angle = 0.0

    @classmethod
    def damage_potential(cls, swing_options):
        """
        This function takes in the basic information about a swing and gives the final resulting model.
        :param swing_options: an object of the SwingOptions class
        :return: swing_result: an object of the SwingResult class
        """
        swing_result = SwingResults()

        ball = Ball(swing_options.ballMass, swing_options.ballRadius)
        bat = Bat(swing_options.batLength, swing_options.batWeight, swing_options.batDiameter, swing_options.batCOR)

        barrel_velocity = [[0.], [swing_options.batSpeed], [0.]]

        # We assume all hits originate from approximately waist height above the front center point of the plate (~3 ft)
        initial_position = np.array([[0.], [0.], [0.91]])

        swung_bat = SwungBat(bat, barrel_velocity, swing_options.rotation_matrix_bat_to_global)

        pitched_ball = PitchedBall.from_incoming_info(ball, swing_options.pitch_speed, swing_options.pitch_angle,
                                                      swing_options.pitch_spin_rate)

        # We first calculate the ideal_velo_flight_results based upon the bat making direct impact with the ball.
        velo_ball_flight_initial_conditions = BallFlightInitialConditions.from_impact(pitched_ball, swung_bat, 0,
                                                                                      initial_position)

        ideal_velo_flight_results = velo_ball_flight_initial_conditions.process(0.01)
        # Note that the above interval of 0.01 was chosen arbitrarily.

        # Establish global frame direction vectors
        global_up = np.array([[0], [0], [1]])
        global_down = np.array([[0], [0], [-1]])
        global_forward = np.array([[0], [1], [0]])

        # Transform those vectors into the bat frame
        bat_up = np.matmul(swung_bat.rotation_matrix_global_to_bat, global_up)
        bat_up[0] = 0  # We remove the x component from each reference vector. We define the x-axis as along the bat,
        #                   so by removing this component, we ensure that each vector is perpendicular to the bat,
        #                   which is how we define impact direction.
        bat_down = np.matmul(swung_bat.rotation_matrix_global_to_bat, global_down)
        bat_down[0] = 0
        bat_forward = np.matmul(swung_bat.rotation_matrix_global_to_bat, global_forward)
        bat_forward[0] = 0

        # Make the transformed vectors into unit vectors
        bat_up = np.multiply(bat_up, 1 / la.norm(bat_up))
        bat_down = np.multiply(bat_down, 1 / la.norm(bat_down))
        bat_forward = np.multiply(bat_forward, 1 / la.norm(bat_forward))

        # Here we take the cross product of the down vector and the forward vector. The resulting vector should point
        #   upwards along the bat (from the end to the handle).
        direction_along_bat = np.cross(np.transpose(bat_down), np.transpose(bat_forward))

        # We assume the direction_along_bat is correct, and points upwards along the bat. However, as a protective
        #   measure, we include this loop. If the vector (which should be [1, 0, 0]) winds up pointing the wrong
        #   direction ([-1, 0, 0]), we want to iterate through the angles in the negative direction so that we don't
        #   get an impact direction that points backwards, towards the backstop.

        direction_correction = 1.0
        if direction_along_bat[0][0] < -0.5:  # We use -0.5 to be safe. The first index should be very close to 1,
            #                                    but just in case something went wrong, we give room for error
            #                                    with -0.5.
            direction_correction = -1.0

        possible_ball_results = []  # We create an array of all the possible results in case we want to play the
        #                                "what-if" game, so that we don't have to reprocess everything, we can simply
        #                                 call a results that occurred at a specific impact angle.
        max_dist = 0.0
        index_max_dist = 0
        rotation_matrix_to_bat_angle = np.zeros([3, 3])

        for i in range(180):
            angle = (direction_correction * i) * DEG_TO_RAD
            offset = ball.radius * np.cos(i * DEG_TO_RAD)  # This gets the magnitude of the vertical distance between
            #                                                    the impact location of the ball and the COM of the bat.

            # Here we create a rotation matrix from the generic bat to the specific angle orientation.
            rotation_matrix_to_bat_angle[0][0] = 1
            rotation_matrix_to_bat_angle[1][1] = np.cos(angle)
            rotation_matrix_to_bat_angle[1][2] = -1.0 * np.sin(angle)
            rotation_matrix_to_bat_angle[2][1] = np.sin(angle)
            rotation_matrix_to_bat_angle[2][2] = np.cos(angle)

            bat_frame_direction_impact = np.matmul(rotation_matrix_to_bat_angle, bat_up)
            swung_bat.bat_frame_impact_direction = bat_frame_direction_impact
            global_impact_direction = np.matmul(swung_bat.rotation_matrix_bat_to_global, bat_frame_direction_impact)
            swung_bat.global_impact_direction = global_impact_direction

            dist_ball_flight_initial_conditions = BallFlightInitialConditions.from_impact(pitched_ball, swung_bat, offset,
                                                                                          initial_position)
            flight_results = dist_ball_flight_initial_conditions.process(0.01)
            # Note that the above interval of 0.01 was chosen arbitrarily.

            possible_ball_results.append(flight_results)
            distance = la.norm(flight_results.flight_position[-1])

            # the horizontal distance traveled is just the magnitude of this vector because the z component will be zero
            #   due to the nature of the flight_model loop (it terminates when z = 0, or when the ball hits the ground)

            if distance > max_dist:
                index_max_dist = i
                max_dist = distance

        ideal_dist_flight_results = possible_ball_results[index_max_dist]

        # ***** Ideal Exit Conditions Calculations *****
        # Note: We call this "ideal" because these calculations at the moment are based upon an ideal iteration. These
        #           are the exit conditions that lead to the greatest total distance traveled, not necessarily the exit
        #           conditions that actually occurred.

        exit_velocity = ideal_velo_flight_results.flight_velocity[0]

        exit_speed = np.linalg.norm(exit_velocity)

        # Launch Angle Calculation
        xy_component = np.linalg.norm([exit_velocity[0], exit_velocity[1]])
        z_component = exit_velocity[2]
        launch_angle = np.arctan(z_component / xy_component)

        # Heading Angle Calculation
        x_component = exit_velocity[0]
        y_component = exit_velocity[1]
        heading_angle = np.arctan2(y_component, x_component)
        # This heading angle is given with 0 being along the x axis, but we need to think of center field as our 0
        #   degree line. The following adjustment does that, where positive angles are towards right field, and negative
        #   angles are towards left.
        heading_angle_field = (np.pi / 2) - heading_angle

        swing_result.flight_position = ideal_dist_flight_results.flight_position
        swing_result.max_distance = max_dist * M_TO_FT
        swing_result.ideal_exit_speed = exit_speed * MS_TO_MPH
        swing_result.ideal_launch_angle = launch_angle * DEG_TO_RAD
        swing_result.ideal_heading_angle = heading_angle_field * DEG_TO_RAD

        # Plotting Functions
        if swing_options.playType == 'SLOW_PITCH_SOFTBALL' or \
                swing_options.playType == 'FAST_PITCH_SOFTBALL':
            softball = True
        else:
            softball = False
        bfp.plot_3d_ball_path([ideal_dist_flight_results], ['Max Dist - ' + str(np.round(max_dist * M_TO_FT, 2)) +
                                                            ' FT'], None, None, False, softball)

        return swing_result


if __name__ == '__main__':
    # pitchoptions = PitchOptions('FAST_PITCH_SOFTBALL', 40., 900., [1, 0, 0], 4.0, 2.0, True)
    # pitchresult = PitchResults.throw_ball(pitchoptions)

    swingoptions = SwingOptions('FAST_PITCH_SOFTBALL', 33.0, 24.0, 55.0, -10.0, 'PITCH')
    swingresult = SwingResults.damage_potential(swingoptions)

