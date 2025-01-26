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
        bob = Dot(point=pendulum.get_end(), color=RED)  # Bob at the end of the pendulum

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
