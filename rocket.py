import numpy as np

# this early in the repo history, this is all artifacts of what I pasted
# from another project

class rocketClass:
    """
    This the the rocket class. It contains all the relevant parameters and
    functions to simulate dynamics and observation.
    """
    def __init__(self):
        # this is the constructor

        # establish constants (see notes to make sense of these)
        # when structured for monte carlo, many of these will be sampled
        self.g      = 9.81      # m/s^2     (earth's gravity)
        self.h_0    = 1188.72   # m         (3900ft, altitude of Las Cruces, NM)
        self.hd_0   = 296       # m/s       (post-burn velocity)                REFINE
        self.rho_0  = 1.036     # kg/m^3    (for 3900ft, 70degF, 15% humidity)  NOTE: could search for a formula for this to be able to vary temp, humidity, etc
        self.T_0    = 294.261   # K         (70degF in Kelvin)                  NOTE: could use formula for this to vary temp
        self.alpha  = 0.0065    # K/m       (low altitude temperature rate)
        self.n      = 5.2561    # unitless  (gas constant)
        self.m      = 20.5      # kg        (post-burn aircraft mass)           REFINE

        self.CD_b   = 0.6       # unitless  (drag base coefficient)             REFINE
        self.CD_s   = 0.02      # unitless  (drag slope CD/angle rate)          FIX FIX FIX (assumes linear relationship)
        self.A      = 0.03      # m^2       (reference area)                    FIX FIX FIX

        # initialize current states
        self.h      = self.h_0
        self.hd     = self.hd_0

        # initialize history vectors
        self.h_all  = np.array([])
        self.hd_all = np.array([])


    def propagateStates(self, dt, theta):
        # propagate full nonlinear equations of motion
        # NOTE: try runge-kutta continuous time to minimize numerical error

        # use states, inputs, and constants to find the derivatives of states
        # first find the drag coefficient based on the current angle theta
        # then use the change in height to find the current air density
        # then find the overall drag and the subsequent acceleration

        CD      = self.CD_b + self.CD_s*theta;

        h_diff  = max(self.h - self.h_0, 0)
        ratio   = (self.T_0 - self.alpha*h_diff)/self.T_0
        rho     = self.rho_0 * ratio**(self.n-1)

        drag    = 0.5*rho*(self.hd**2)*CD*self.A
        hdd     = -self.g - (drag/self.m)*np.sign(self.hd)

        # update the states with crude foward difference
        # could be improved with numerical method
        self.hd = self.hd + dt*hdd
        self.h  = self.h  + dt*self.hd


        # note, the append gets slow, figure out good increment system and preallocate
        self.h_all  = np.append(self.h_all, np.array(self.h))
        self.hd_all = np.append(self.hd_all, np.array(self.hd))
