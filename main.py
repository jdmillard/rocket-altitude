import numpy as np
import time

from plotter import Plotter


h = 0
h_list = np.array([])

v_up = 30
decel = 10

dt = 0.01
start_time = 0
final_time = 1
times = np.arange(start_time,final_time+dt,dt)
t_list = np.array([])

plotter = Plotter(final_time=final_time, plot_real_time=True)

for t in times:

    h = h + v_up*dt
    h_list = np.append(h_list, np.array(h))
    t_list = np.append(t_list, np.array(t))

    v_up = v_up - decel*dt

    plotter.updateItems(h_list, t_list, t, time.time())








print('end')
