# import core libraries
import numpy as np
import time

# import custom classes
from plotter import livePlotter
from rocket import rocketClass

# import signal handler (for clean ctrl+C exiting)
import signal
import sys

# define operations during shutdown
def signal_handler(signal, frame):
    print('\nterminating simulation')
    sys.exit(0)

# setup signal handler specifying callback function
signal.signal(signal.SIGINT, signal_handler)



# set dt, start, end, and create array of all times
dt = 0.01
start_time = 0
final_time = 20
times = np.arange(start_time+dt,final_time+dt,dt)

# instantiate plotter and rocket classes
plotter = livePlotter(final_time=final_time, plot_real_time=True)
rocket = rocketClass(times)

# perform simulation
for t in times:

    rocket.propagateStates(dt, 0)
    plotter.updateItems(rocket, t, time.time())





print('end')
