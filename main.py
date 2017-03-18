import numpy as np
import time

from plotter import livePlotter
from rocket import rocketClass


# initialize history vectosr for altitude, velocity, and time
h_list = np.array([])
hd_list = np.array([])
t_list = np.array([])

# set dt, start, end, and create array of all times
dt = 0.01
start_time = 0
final_time = 20
times = np.arange(start_time,final_time+dt,dt)

# instantiate plotter and rocket classes
plotter = livePlotter(final_time=final_time, plot_real_time=True)
rocket = rocketClass()

# perform simulation
for t in times:

    rocket.propagateStates(dt, 0)

    h = rocket.h
    hd = rocket.hd

    h_list  = np.append(h_list, np.array(h))
    hd_list = np.append(h_list, np.array(h))
    t_list  = np.append(t_list, np.array(t))

    plotter.updateItems(h_list, t_list, t, time.time())





print('end')
