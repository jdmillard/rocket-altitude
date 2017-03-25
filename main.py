# import core libraries
import numpy as np
import time

# import custom classes
from rocket import rocketClass
from plotter import livePlotter
from lcm import lcm_mult as lcm

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

'''

# set the rates of each piece of simulation
hz_sim = 100    # (Hz) sim rate
hz_fil = 50     # (Hz) filter and control rate
hz_plt = 60     # (Hz) plot rate

# get the rate that represents the least common multiple
hz_lcm = lcm([hz_sim, hz_fil, hz_plt])

# set dt, start, end, and create array of all times
dt = 0.01
start_time = 0
final_time = 20
times = np.arange(start_time+dt,final_time+dt,dt)

# instantiate plotter and rocket classes
rocket  = rocketClass(times)
plotter = livePlotter(rocket, final_time=final_time, plot_real_time=True)

# perform simulation
for t in times:

    rocket.propagateStates(dt, 0)
    plotter.updateItems(rocket, t, time.time())





print('end')
