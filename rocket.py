import numpy as np


class rocketClass:
    """
    This the the rocket class. It contains all the relevant parameters and
    functions to simulate dynamics and observation.
    """
    def __init__(self, times):
        # this is the constructor

        # establish constants (see notes to make sense of these)
        # when structured for monte carlo, many of these will be sampled
        # best to consider these beyond our control
        self.g      = 9.81      # m/s^2     (earth's gravity)
        self.h_0    = 1188.72   # m         (3900ft, Las Cruces, NM)            NOTE: we need to take into account the amount of altitude gain during burn
        self.h_f    = 3048      # m         (10000ft, final destination)
        self.hd_0   = 296       # m/s       (post-burn velocity)                REFINE
        self.rho_0  = 1.036     # kg/m^3    (for 3900ft, 70degF, 15% humidity)  NOTE: could search for a formula for this to be able to vary temp, humidity, etc
        self.T_0    = 294.261   # K         (70degF in Kelvin)                  NOTE: could use formula for this to vary temp
        self.alpha  = 0.0065    # K/m       (low altitude temperature rate)
        self.n      = 5.2561    # unitless  (gas constant)
        self.m      = 20.5      # kg        (post-burn aircraft mass)           REFINE

        self.CD_b   = 0.6       # unitless  (drag base coefficient)             REFINE
        self.CD_s   = 0.02      # unitless  (drag slope CD/angle rate)          FIX FIX FIX (assumes linear relationship)
        self.A      = 0.025     # m^2       (reference area)                    FIX FIX FIX

        # initialize current states
        self.h      = self.h_0
        self.hd     = self.hd_0

        # initialize history vectors
        self.h_all  = np.empty(times.size+1)
        self.hd_all = np.empty(times.size+1)
        self.t_all  = np.empty(times.size+1)

        self.h_all[0]   = self.h
        self.hd_all[0]  = self.hd
        self.t_all[0]   = 0
        self.i          = 1 # index for next iteration

        # establish the reference trajectory
        CD_ref = self.CD_b * 1.2; # NOTE: in the future we won't know CD_b, this is only a temporary assignment
        self.refTrajectory(CD_ref)


    def propagateStates(self, dt, theta):
        # propagate full nonlinear equations of motion
        # NOTE: try runge-kutta continuous time to minimize numerical error

        # use states, inputs, and constants to find the derivatives of states
        # first find the drag coefficient based on the current angle theta
        # then use the change in height to find the current air density
        # then find the overall drag and the subsequent acceleration

        CD      = self.CD_b + self.CD_s*theta;

        h_diff  = self.h - self.h_0
        ratio   = (self.T_0 - self.alpha*h_diff)/self.T_0
        rho     = self.rho_0 * ratio**(self.n-1)

        drag    = 0.5*rho*(self.hd**2)*CD*self.A
        hdd     = -self.g - (drag/self.m)*np.sign(self.hd)

        # update the states with crude foward difference
        # could be improved with numerical method
        self.hd = self.hd + dt*hdd
        self.h  = self.h  + dt*self.hd

        # update the history vectors
        self.h_all[self.i]  = self.h
        self.hd_all[self.i] = self.hd
        self.t_all[self.i]  = self.t_all[self.i-1]+dt

        # increment index
        self.i = self.i + 1


    def refTrajectory(self, CD_ref):
        # use pessimistic parameters to define a reference trajectory (hd vs h)
        # this is done starting from sea level in order to be more complete
        # vary hd_0_ref until the open loop trajectory arrives at desired h

        # in the future, if we have a good on-line drage estimation, this
        # reference trajectory could be updated during flight, but for now
        # we'll use pessimistic settings

        hd_0_ref    = self.hd_0 # guess for starting velocity at sea level
        dt          = 0.01      # only for reference trajectory generation
        disparity   = 100       # arbitrarily high start for while loop
        thresh      = 0.01      # threshold for desired altitude

        while abs(disparity) > thresh:
            # so long as the disparity is too high, keep guessing
            # reset h, hd, and their histories
            h       = 0
            hd      = hd_0_ref
            h_ref   = np.empty([0])
            hd_ref  = np.empty([0])

            # run the trajectory with the current parameters
            while hd > 0:
                # find air density at this altitude
                h_diff  = h - self.h_0
                ratio   = (self.T_0 - self.alpha*h_diff)/self.T_0
                rho     = self.rho_0 * ratio**(self.n-1)

                # get drag and deceleration
                drag    = 0.5*rho*(hd**2)*CD_ref*self.A
                hdd     = -self.g - (drag/self.m)*np.sign(hd)

                # update the states with crude foward difference
                # could be improved with numerical method
                hd = hd + dt*hdd
                h  = h  + dt*hd

                # update the history vectors
                # the if check is to avoid writing early, saving time
                if disparity < thresh*10:
                    h_ref  = np.append(h_ref,  [h] )
                    hd_ref = np.append(hd_ref, [hd])

            # update the inital hd based on arrival trajectory
            disparity = self.h_f - h
            hd_0_ref = hd_0_ref + 0.1*(disparity)

        # save the final reference trajectory as a class member
        self.h_ref  = h_ref
        self.hd_ref = hd_ref
