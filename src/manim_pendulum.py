from manim import *

class OscillatingPendulum(Scene):
    def construct(self):
        # Constants
        L = 2  # Length of the pendulum
        g = 9.81  # Acceleration due to gravity
        theta0 = 30 * DEGREES  # Initial angle (30 degrees)
        omega = 0  # Initial angular velocity
        dt = 0.05  # Time step
        t_max = 4   # Total simulation time

        # Set background to transparent
        self.camera.background_color = "#fdf6e3"

        # Ceiling
        ceiling = Line(LEFT * 3, RIGHT * 3, color="#586e75").shift(UP * 0.5)
        self.add(ceiling)

        # Pivot point (attached to the ceiling)
        pivot = ceiling.get_center()# + UP * 0.5

        # Pendulum setup
        pendulum = Line(pivot, pivot + L * DOWN, color="#586e75")  # Pendulum rod with custom color
        bob = Dot(point=pendulum.get_end(), color=RED, radius = 0.2)  # Bob at the end of the pendulum

        # Angle label
        angle_label = always_redraw(lambda: MathTex(
            r"\theta = {:.2f}^\circ".format(np.degrees(theta0)),
            color="#586e75"
        ).next_to(pivot, UP, buff=0.5))

        # Add pendulum and angle label to the scene
        self.add(pendulum, bob, angle_label)

        # Function to update the pendulum's position
        def update_pendulum(pendulum, bob, dt):
            nonlocal theta0, omega
            alpha = -(g / L) * np.sin(theta0)  # Angular acceleration
            omega += alpha * dt  # Update angular velocity
            theta0 += omega * dt  # Update angle
            # Update pendulum and bob positions
            new_end = pivot + L * np.sin(theta0) * RIGHT + L * np.cos(theta0) * DOWN
            pendulum.put_start_and_end_on(pivot, new_end)
            bob.move_to(new_end)

        # Animate the pendulum in a loop
        for _ in range(4):  # Loop 3 times (you can adjust this or remove the loop for infinite animation)
            self.play(
                UpdateFromFunc(
                    VGroup(pendulum, bob),
                    lambda mob: update_pendulum(pendulum, bob, dt)
                ),
                rate_func=linear,
                run_time=t_max
            )
            # Reset the pendulum to its initial state for the next loop
            theta0 = 30 * DEGREES
            omega = 0
            new_end = pivot + L * np.sin(theta0) * RIGHT + L * np.cos(theta0) * DOWN
            pendulum.put_start_and_end_on(pivot, new_end)
            bob.move_to(new_end)

        # Hold the final frame
        self.wait(2)


class DoublePendulumComparison(Scene):
    def construct(self):
        # Constants
        L = 2  # Length of the pendulum
        g = 9.81  # Acceleration due to gravity
        theta0_1 = 170 * DEGREES  # Initial angle for pendulum 1 (exact)
        theta0_2 = 170 * DEGREES  # Initial angle for pendulum 2 (approximate)
        omega1 = 0  # Initial angular velocity for pendulum 1
        omega2 = 0  # Initial angular velocity for pendulum 2
        dt = 0.05  # Time step
        t_max = 4   # Total simulation time

        # Set background to transparent
        self.camera.background_color = "#fdf6e3"

        # Ceiling
        ceiling = Line(LEFT * 3, RIGHT * 3, color="#586e75").shift(UP * 0.5)
        self.add(ceiling)

        # Pivot points (attached to the ceiling)
        pivot1 = ceiling.get_center() + LEFT * 1.8
        pivot2 = ceiling.get_center() + RIGHT * 1.8

        # Pendulum 1 setup (exact)
        pendulum1 = Line(pivot1, pivot1 + L * DOWN, color="#586e75")  # Pendulum rod with custom color
        bob1 = Dot(point=pendulum1.get_end(), color=RED, radius=0.2)  # Bob at the end of the pendulum

        # Pendulum 2 setup (approximate)
        pendulum2 = Line(pivot2, pivot2 + L * DOWN, color="#586e75")  # Pendulum rod with custom color
        bob2 = Dot(point=pendulum2.get_end(), color=BLUE, radius=0.2)  # Bob at the end of the pendulum

        # Angle labels
        angle_label1 = always_redraw(lambda: MathTex(
            r"\theta_1 = {:.2f}^\circ".format(np.degrees(theta0_1)),
            color="#586e75"
        ).next_to(pivot1, UP, buff=0.5))

        angle_label2 = always_redraw(lambda: MathTex(
            r"\theta_2 = {:.2f}^\circ".format(np.degrees(theta0_2)),
            color="#586e75"
        ).next_to(pivot2, UP, buff=0.5))

        # Add pendula and angle labels to the scene
        self.add(pendulum1, bob1, angle_label1, pendulum2, bob2, angle_label2)

        # Function to update pendulum 1's position (exact formula)
        def update_pendulum1(pendulum, bob, dt):
            nonlocal theta0_1, omega1
            alpha = -(g / L) * np.sin(theta0_1)  # Exact angular acceleration
            omega1 += alpha * dt  # Update angular velocity
            theta0_1 += omega1 * dt  # Update angle
            # Update pendulum and bob positions
            new_end = pivot1 + L * np.sin(theta0_1) * RIGHT + L * np.cos(theta0_1) * DOWN
            pendulum.put_start_and_end_on(pivot1, new_end)
            bob.move_to(new_end)

        # Function to update pendulum 2's position (approximate formula)
        def update_pendulum2(pendulum, bob, dt):
            nonlocal theta0_2, omega2
            alpha = -(g / L) * theta0_2  # Small-angle approximation for angular acceleration
            omega2 += alpha * dt  # Update angular velocity
            theta0_2 += omega2 * dt  # Update angle
            # Update pendulum and bob positions
            new_end = pivot2 + L * np.sin(theta0_2) * RIGHT + L * np.cos(theta0_2) * DOWN
            pendulum.put_start_and_end_on(pivot2, new_end)
            bob.move_to(new_end)

        # Animate the pendula in a loop
        for _ in range(2):  # Loop 2 times (you can adjust this or remove the loop for infinite animation)
            self.play(
                UpdateFromFunc(
                    VGroup(pendulum1, bob1),
                    lambda mob: update_pendulum1(pendulum1, bob1, dt)
                ),
                UpdateFromFunc(
                    VGroup(pendulum2, bob2),
                    lambda mob: update_pendulum2(pendulum2, bob2, dt)
                ),
                rate_func=linear,
                run_time=t_max
            )
            # Reset the pendula to their initial states for the next loop
            theta0_1 = 170 * DEGREES
            theta0_2 = 170 * DEGREES
            omega1 = 0
            omega2 = 0
            new_end1 = pivot1 + L * np.sin(theta0_1) * RIGHT + L * np.cos(theta0_1) * DOWN
            new_end2 = pivot2 + L * np.sin(theta0_2) * RIGHT + L * np.cos(theta0_2) * DOWN
            pendulum1.put_start_and_end_on(pivot1, new_end1)
            pendulum2.put_start_and_end_on(pivot2, new_end2)
            bob1.move_to(new_end1)
            bob2.move_to(new_end2)

        # Hold the final frame
        self.wait(2)
