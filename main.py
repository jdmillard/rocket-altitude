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

if (True):
    # set dt, start, end, and create array of all times
    dt = 0.01
    start_time = 0
    final_time = 15
    times = np.arange(start_time+dt,final_time+dt,dt)

    # instantiate plotter and rocket classes
    rocket  = rocketClass(times, 1)
    plotter = livePlotter(rocket, final_time=final_time, plot_real_time=False)

    # perform simulation
    for t in times:

        rocket.setControl(dt)
        rocket.setControl2(dt)
        rocket.propagateStates(dt)
        plotter.updateItems(rocket, t, time.time())
        # need to implement adaptive plotting framerate
        # because updating the data in the plotting objects takes a long time
        # once the vectors are large: this is why it drops below real-time





def oneRun(Qin, hh):
    # set dt, start, end, and create array of all times
    dt = 0.01
    start_time = 0
    final_time = 6
    times = np.arange(start_time+dt,final_time+dt,dt)

    # instantiate plotter and rocket classes
    rocket  = rocketClass(times, hh)
    rocket.Q = Qin
    #plotter = livePlotter(rocket, final_time=final_time, plot_real_time=False)

    # perform simulation
    for t in times:

        rocket.setControl(dt)
        rocket.setControl2(dt)
        rocket.propagateStates(dt)
        #plotter.updateItems(rocket, t, time.time())
        # need to implement adaptive plotting framerate
        # because updating the data in the plotting objects takes a long time
        # once the vectors are large: this is why it drops below real-time

    return rocket.cum_error





Q = np.array([[0.4659 , 0.001 , 0.001 , 0.001 , 0.004 , 0.001 , 0.001 ],
              [0.001  , 15863 , 0.001 , 0.001 , 0.000125 , 0.001 , 0.001 ],
              [0.001  , 0.001 , 0.00036, 0.001, 0.016 , 0.001 , 0.001 ],
              [0.001  , 0.001 , 0.001 , 0.21  , 0.0005 , 0.001 , 0.001 ],
              [0.004  , 0.000125 , 0.016 , 0.0005 , 0.00336, 0.001 , 0.001 ],
              [0.001  , 0.001 , 0.001 , 0.001 , 0.001 , 0.001 , 0.001 ],
              [0.001  , 0.001 , 0.001 , 0.001 , 0.001 , 0.001 , 0.001 ]])

scale = 1.05


# currently trying to converge on the CD_bar with all cross terms
# then it would be nice to open it up to all terms and run overnight

# eventually, make the error the sum of all states and cycle through all
# Q elements, using a larger amount of runs for a better average
# (maybe)

# put the tuning logic in a separate .py
# idea - save Q to file and load it for improved management?
# have each run randomly sample a few parameters

# need logic for CD convergence since it becomes unobservable at low speeds
# perhaps bounds on the estimation as well to prevent filter breakage
# wait until after subsonic speeds are reached, capture before speed drops too near-unobservable levels, lpf
# try catch on F = linalg.expm(self.A*dt)

'''
# latest overnight run, not looked at yet:
[[  4.65900000e-01   1.00000000e-03   1.00000000e-03   1.00000000e-03    4.20000000e-03]
 [  1.00000000e-03   1.58630000e+04   1.00000000e-03   1.00000000e-03    1.44703125e-04]
 [  1.00000000e-03   1.00000000e-03   3.60000000e-04   1.00000000e-03    1.45124717e-02]
 [  1.00000000e-03   1.00000000e-03   1.00000000e-03   2.10000000e-01    5.78812500e-04]
 [  4.20000000e-03   1.44703125e-04   1.45124717e-02   5.78812500e-04    4.08410100e-03]]
'''
batches = 50
iii = 0
while iii < batches:
    iii = iii + 1
    if (iii>20):
        scale = 1.01
    h_list = np.arange(4,5) # states to consider
    for h in h_list:
        print('-------')
        print("state", h+1)

        i_list = np.arange(0,5) # cross terms to consider
        for i in i_list:
            # h is the current state
            # i is the cross term

            Q_dec = Q.copy()
            if (i==h):
                Q_dec[h,i] = Q_dec[h,i] * (1/scale)
            else:
                Q_dec[h,i] = Q_dec[h,i] * (1/scale)
                Q_dec[i,h] = Q_dec[i,h] * (1/scale)
            err1 = oneRun(Q_dec, h)
            err2 = oneRun(Q_dec, h)
            err3 = oneRun(Q_dec, h)
            err4 = oneRun(Q_dec, h)
            err5 = oneRun(Q_dec, h)
            err_dec = (err1 + err2 + err3 + err4 + err5)/5

            Q_mid = Q.copy()
            err1 = oneRun(Q_mid, h)
            err2 = oneRun(Q_mid, h)
            err3 = oneRun(Q_mid, h)
            err4 = oneRun(Q_mid, h)
            err5 = oneRun(Q_mid, h)
            err_mid = (err1 + err2 + err3 + err4 + err5)/5


            Q_inc = Q.copy()
            if (i==h):
                Q_inc[h,i] = Q_inc[h,i] * scale
            else:
                Q_inc[h,i] = Q_inc[h,i] * scale
                Q_inc[i,h] = Q_inc[i,h] * scale
            err1 = oneRun(Q_inc, h)
            err2 = oneRun(Q_inc, h)
            err3 = oneRun(Q_inc, h)
            err4 = oneRun(Q_inc, h)
            err5 = oneRun(Q_inc, h)
            err_inc = (err1 + err2 + err3 + err4 + err5)/5

            print(err_mid)
            if ((err_inc < err_dec) and (err_inc < err_mid)):
                Q = Q_inc.copy()
                print(err_inc)
                #print(Q)
            elif ((err_dec < err_inc) and (err_dec < err_mid)):
                Q = Q_dec.copy()
                print(err_dec)
                #print(Q)
            else:
                print(err_mid)

            print(Q[0:5,0:5])


print(Q[0:5,0:5])
