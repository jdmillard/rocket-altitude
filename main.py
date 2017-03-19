import numpy as np
import time

from plotter import livePlotter
from rocket import rocketClass


# initialize history vector for time
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

    t_list  = np.append(t_list, np.array(t))

    rocket.propagateStates(dt, 0)

    plotter.updateItems(rocket, t_list, t, time.time())





print('end')
