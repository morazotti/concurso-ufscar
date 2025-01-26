import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from scipy.interpolate import interp1d

transparency = True
dpi = 300
color = "#586e75"


def set_ticks(plot):
    plot.xticks(
        [i * 0.5 * np.pi for i in range(9)],
        [
            "0",
            r"$T/4$",
            r"$T/2$",
            r"$3T/4$",
            r"$T$",
            r"$5T/4$",
            r"$3T/2$",
            r"$7T/4$",
            r"$2T$",
        ],
    )
    plot.xlabel(r"$t$")
    plot.ylabel(r"$\theta$ (rad)")
    # set frame color
    plot.gca().spines["top"].set_color("none")
    plot.gca().spines["right"].set_color("none")
    plot.gca().spines["left"].set_color(f"{color}")
    plot.gca().spines["bottom"].set_color(f"{color}")
    # set ticks color
    plot.gca().xaxis.label.set_color(f"{color}")
    plot.gca().yaxis.label.set_color(f"{color}")
    plot.tick_params(axis="x", colors=f"{color}")
    plot.tick_params(axis="y", colors=f"{color}")


class Pendulum:

    def __init__(self, L, theta0, omega0, exact):
        self.L = L
        self.theta0 = theta0
        self.omega0 = omega0
        self.exact = exact
        self.w_approx = np.sqrt(9.81 / self.L)
        self.T = 2 * np.pi / self.w_approx

    def ode_system(self, t, y):
        f, fprime = y
        if self.exact:
            return [fprime, -9.81 * np.sin(f) / self.L]
        return [fprime, -9.81 * f / self.L]

    def solve_ode(self):
        y0 = [self.theta0, self.omega0]

        # Solve the ODE between t = 0 and t = 4π
        t_span = (0, 4 * np.pi / self.w_approx)
        t_eval = np.linspace(t_span[0], t_span[1], 1000)  # Fine grid for interpolation
        solution = solve_ivp(self.ode_system, t_span, y0, t_eval=t_eval, method="RK45")

        # Interpolate the solution
        t_values = solution.t
        f_values = solution.y[0]
        func = interp1d(t_values, f_values, kind="cubic", fill_value="extrapolate")

        return func


L = 10
theta0, omega0 = np.pi / 18, 0
nontheta = Pendulum(L, theta0, omega0, True).solve_ode()
theta = Pendulum(L, theta0, omega0, False)
w = theta.w_approx
theta = theta.solve_ode()
t = np.linspace(0, 4 * np.pi / w, 100)

fig, ax = plt.subplots(dpi=dpi)
set_ticks(plt)
# ax.plot(t, A*np.sin(w*t + phi), color='steelblue', label=r'$\varphi = $'+f'{phi:.2f}')
ax.plot(t, nontheta(t), color="palevioletred", label="Solução exata")
ax.plot(t, theta(t), color="steelblue", label="Solução aproximada", linestyle="--")
ax.legend(loc="upper right", fontsize=8)
ax.set_title(
    r"Pêndulo não-linear, $\theta_0 = {:.2f}^\circ, \omega_0 = {:.2f}$".format(
        180 * theta0 / np.pi, omega0
    )
)
plt.savefig("../img/nonlinear_pendulum_close.png", transparent=transparency)

theta0, omega0 = np.pi / 3, 0
nontheta = Pendulum(L, theta0, omega0, True).solve_ode()
theta = Pendulum(L, theta0, omega0, False)
w = theta.w_approx
theta = theta.solve_ode()
t = np.linspace(0, 4 * np.pi / w, 100)

fig, ax = plt.subplots(dpi=dpi)
set_ticks(plt)
# ax.plot(t, A*np.sin(w*t + phi), color='steelblue', label=r'$\varphi = $'+f'{phi:.2f}')
ax.plot(t, nontheta(t), color="palevioletred", label="Solução exata")
ax.plot(t, theta(t), color="steelblue", label="Solução aproximada", linestyle="--")
ax.legend(loc="upper right", fontsize=8)
ax.set_title(
    r"Pêndulo não-linear, $\theta_0 = {:.2f}^\circ, \omega_0 = {:.2f}$".format(
        180 * theta0 / np.pi, omega0
    )
)
plt.savefig("../img/nonlinear_pendulum_far.png", transparent=transparency)

theta0, omega0 = 0.99*np.pi, 0
nontheta = Pendulum(L, theta0, omega0, True).solve_ode()
theta = Pendulum(L, theta0, omega0, False)
w = theta.w_approx
theta = theta.solve_ode()
t = np.linspace(0, 4 * np.pi / w, 100)

fig, ax = plt.subplots(dpi=dpi)
set_ticks(plt)
# ax.plot(t, A*np.sin(w*t + phi), color='steelblue', label=r'$\varphi = $'+f'{phi:.2f}')
ax.plot(t, nontheta(t), color="palevioletred", label="Solução exata")
ax.plot(t, theta(t), color="steelblue", label="Solução aproximada", linestyle="--")
ax.legend(loc="upper right", fontsize=8)
ax.set_title(
    r"Pêndulo não-linear, $\theta_0 = {:.2f}^\circ, \omega_0 = {:.2f}$".format(
        180 * theta0 / np.pi, omega0
    )
)
plt.savefig("../img/nonlinear_pendulum_upsidedown.png", transparent=transparency)
