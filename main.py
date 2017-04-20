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
    final_time = 9
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





Q = np.array([[0.35   , 0.001 , 0.001 , 0.001 , 0.001 , 0.001 , 0.01  ],
              [0.001  , 1100  , 0.001 , 0.001 , 0.001 , 0.001 , 0.01  ],
              [0.001  , 0.001 , 0.1   , 0.001 , 0.001 , 0.001 , 0.01  ],
              [0.001  , 0.001 , 0.001 , 20    , 0.001 , 0.001 , 0.01  ],
              [0.001  , 0.001 , 0.001 , 0.001 , 0.001 , 0.001 , 0.01  ],
              [0.001  , 0.001 , 0.001 , 0.001 , 0.001 , 0.001 , 0.01  ],
              [0.001  , 0.001 , 0.001 , 0.001 , 0.001 , 0.001 , 0.01  ]])

scale = 1.1


# added lpf to hd, now tune for least squares
# get h and h_dot tuned up as best as possible

# then open up theta and get these 3 states tuned up nicely
# then open up theta_d and get theta+theta_d tuned nicely together
# then get these 4 states tuned nicely together

h_list = np.arange(0,1) # states to consider
for h in h_list:
    print('-------')
    print("state", h+1)

    i_list = np.arange(0,1) # cross terms to consider
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


print(Q[0:4,0:4])
