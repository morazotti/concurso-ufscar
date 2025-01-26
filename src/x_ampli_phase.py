import matplotlib.pyplot as plt
import numpy as np


transparency = True
dpi = 300
color = "#586e75"
w, A, phi = 1, 1, 0
v0, x0 = -w*A*np.sin(phi), A*np.cos(phi)
t = np.linspace(0, 4*np.pi/w, 100)


def set_ticks(plot, energy = False):
    plot.xticks([i*0.5*np.pi for i in range(9)], ['0', r'$T/4$', r'$T/2$', r'$3T/4$', r'$T$', r'$5T/4$', r'$3T/2$', r'$7T/4$', r'$2T$'])
    plot.xlabel(r'$t$')
    if not energy:
        plot.yticks([-1, 0, 1], ['-A', '0', 'A'])
        plot.ylabel(r'$x(t)$')
    # set frame color
    plot.gca().spines['top'].set_color('none')
    plot.gca().spines['right'].set_color('none')
    plot.gca().spines['left'].set_color(f'{color}')
    plot.gca().spines['bottom'].set_color(f'{color}')
    # set ticks color
    plot.gca().xaxis.label.set_color(f'{color}')
    plot.gca().yaxis.label.set_color(f'{color}')
    plot.tick_params(axis='x', colors=f'{color}')
    plot.tick_params(axis='y', colors=f'{color}')


fig, ax = plt.subplots(dpi=dpi)
set_ticks(plt)
ax.plot(t, A*np.sin(w*t + phi), color='steelblue')
plt.savefig('../img/x_ampli_phase.png', transparent=transparency)


fig, ax = plt.subplots(dpi=dpi)
set_ticks(plt)
ax.plot(t, A*np.sin(w*t + phi), color='steelblue', label=r'$\varphi = $'+f'{phi:.2f}')
ax.plot(t, A*np.sin(w*t + phi + np.pi/2.5), color='palevioletred', label=r'$\varphi = $' + f'{phi + np.pi/2.5:.2f}')
ax.legend(loc='upper right', fontsize=8)
plt.savefig('../img/x_change_phase.png', transparent=transparency)


fig, ax = plt.subplots(dpi=dpi)
set_ticks(plt)
ax.plot(t, A*np.sin(w*t + phi), color='steelblue', label=f'$A = {A}$')
ax.plot(t, A*np.sin(w*t + phi)/1.5, color='palevioletred', label=f'$A = {A/1.5:.2f}$')
ax.legend(loc='upper right', fontsize=8)
plt.savefig('../img/x_change_ampli.png', transparent=transparency)


fig, ax = plt.subplots(dpi=dpi)
set_ticks(plt)
ax.plot(t, x0*np.cos(w*t) + v0*np.sin(w*t), color='steelblue', label=r'$x_0 = $'+f'{x0}, '+r'$v_0 = $'+f'{v0}')
ax.plot(t, (x0 - 0.3)*np.cos(w*t) + (v0 + 1.5)*np.sin(w*t), color='palevioletred', label=r'$x_0 = $'+f'{x0-0.3}, ' +r'$v_0 = $'+f'{v0+1.5}')
ax.legend(loc='upper right', fontsize=8)
plt.savefig('../img/x_change_init_cond.png', transparent=transparency)


fig, ax = plt.subplots(dpi=dpi)
set_ticks(plt)
ax.plot(t, A*np.sin(w*t + phi), color='steelblue', label=r'$\omega = $'+f'{w}')
ax.plot(t, A*np.sin(0.5*w*t + phi), color='palevioletred', label=r'$\omega = $'+f'{0.5*w}')
ax.legend(loc='upper right', fontsize=8)
plt.savefig('../img/x_change_freq.png', transparent=transparency)

# energy stuff
x = A*np.cos(w*t + phi)
v = -A*w*np.sin(w*t + phi)
T = v * v / 2
V = w * w * x * x / 2
E = T + V

fig, ax = plt.subplots(dpi=dpi)
set_ticks(plt, True)
plt.yticks([0, 0.5*E[0], E[0]], ['0', r'$\frac{m\omega^2A^2}{4}$', r'$\frac{m\omega^2A^2}{2}$'])
ax.set_title(r'$x_0 = $'+f'{x0}, '+r'$v_0 = $'+f'{v0}')
ax.plot(t, V, color='steelblue', label='Energia potencial')
ax.plot(t, T, color='palevioletred', label='Energia cinética')
ax.plot(t, E, color='seagreen', label='Energia total')
ax.legend(loc='upper right', fontsize=8)
plt.savefig('../img/energy.png', transparent=transparency)
