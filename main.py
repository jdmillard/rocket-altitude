# import core libraries
import numpy as np
import time

# import custom classes
from rocket import rocketClass
from plotter import livePlotter

# import signal handler (for clean ctrl+C exiting)
import signal
import sys

# define operations during shutdown
def signal_handler(signal, frame):
    print('\nterminating simulation')
    sys.exit(0)

# setup signal handler specifying callback function
signal.signal(signal.SIGINT, signal_handler)


'''
This is main.py which initializes the important objects and governs the timing
of each simulation piece.

The simulation starts immediately after burn has finished.
'''


# set dt, start, end, and create array of all times
dt = 0.01
start_time = 0
final_time = 16
times = np.arange(start_time+dt,final_time+dt,dt)

# instantiate plotter and rocket classes
rocket  = rocketClass(times)
plotter = livePlotter(rocket, final_time=final_time, plot_real_time=False)

# perform simulation
for t in times:

    rocket.setControl(dt)
    rocket.propagateStates(dt)
    plotter.updateItems(rocket, t, time.time())
    # need to implement adaptive plotting framerate
    # because updating the data in the plotting objects takes a long time
    # once the vectors are large: this is why it drops below real-time





print('end')
