import numpy as np
import scipy.linalg as linalg


class rocketClass:
    """
    This the the rocket class. It contains all the relevant parameters and
    functions to simulate dynamics and observation.
    """
    def __init__(self, times):
        # this is the constructor, it runs once upon class instantiation
        # "self." means the variable is a member of this class

        # establish constants
        # when structured for monte carlo, many of these will be sampled
        # best to consider these beyond our control
        self.g      = 9.81      # m/s^2     (earth's gravity)
        self.h_0    = 1188.72   # m         (3900ft, Las Cruces, NM)
        self.h_b    = 893.9     # m         (altitude gained during burn)       REFINE, LARGE VARIANCE
        self.h_f    = 3048      # m         (10000ft, final destination)
        self.hd_0   = 303.4     # m/s       (post-burn velocity)                REFINE, LARGE VARIANCE
        self.rho_0  = 1.036     # kg/m^3    (for 3900ft, 70degF, 15% humidity)  NOTE: could search for a formula for this to be able to vary temp, humidity, etc
        self.T_0    = 294.261   # K         (70degF in Kelvin)                  NOTE: could use formula for this to vary temp
        self.alpha  = 0.0065    # K/m       (low altitude temperature rate)
        self.n      = 5.2561    # unitless  (gas constant)
        self.m      = 18.03     # kg        (post-burn aircraft mass)           REFINE, MEDIUM VARIANCE

        self.CD_b   = 2.1       # unitless  (drag base coefficient)             REFINE 0.6 questionably from openrocket
        self.CD_s   = 0.21      # unitless  (drag slope CD/angle rate)          FIX FIX FIX (arbitrarily GUESSED and assumes linear relationship, needs complete overhaul)
        self.A_ref  = 0.0192    # m^2       (reference area, cross section)
        self.th_max = 60        # deg       (maximum air brake angle = 90deg)   COULD EVENTUALLY CHANGE
        self.th_r   = 45        # deg/s     (air brake rate of change)          0 to 90deg in 2 seconds

        # initialize current states and inputs
        self.h      = self.h_0 + self.h_b
        self.hd     = self.hd_0
        self.th     = 0 # theta dot
        self.thd    = 0 # theta dot
        self.th_cmd = 0
        self.u      = 0

        # initialize state truth and estimate vectors
        self.x_tru  = np.zeros((7,1))
        self.x_hat  = np.zeros((7,1))

        # populate the state truth vector
        self.x_tru[0,0] = self.h   # altitude truth
        self.x_tru[1,0] = self.hd  # altitude change truth
        self.x_tru[2,0] = self.th  # theta truth
        self.x_tru[3,0] = self.thd # theta change truth
        self.x_tru[4,0] = self.CD_b*self.A_ref
        self.x_tru[5,0] = self.CD_s
        self.x_tru[6,0] = self.g

        # initialize the filter matrices
        self.P = np.zeros((7,7)) # covariance matrix, P
        self.A = np.zeros((7,7)) # homgeneous system
        self.A[0,1] = 1     # remaining terms are dynamic
        self.A[1,6] = -1    # remaining terms are dynamic
        self.A[2,3] = 1     # remaining terms are dynamic


        # initialize history vectors for plotting
        self.h_all      = np.empty(times.size+1)    # history of h
        self.hd_all     = np.empty(times.size+1)    # history of h_dot
        self.t_all      = np.empty(times.size+1)    # history of time
        self.e_hd       = np.empty(times.size+1)    # history of ref error
        self.th_all     = np.empty(times.size+1)    # history of theta
        self.th_cmd_all = np.empty(times.size+1)    # history of theta desired
        self.x_tru_all  = np.empty((7,times.size+1))
        self.x_hat_all  = np.empty((7,times.size+1))

        # fill first element of history vectors
        self.h_all[0]       = self.h
        self.hd_all[0]      = self.hd
        self.th_all[0]      = self.th
        self.th_cmd_all[0]  = self.th_cmd
        self.t_all[0]       = 0
        self.x_tru_all[:,0] = self.x_tru.ravel()

        # index for next iteration
        self.i              = 1

        # establish the reference trajectory
        CD_ref = self.CD_b * 1.2                                                # NOTE: in the future we won't know CD_b, this is only a temporary assignment
        self.refTrajectory(CD_ref)
        self.integrator = 0;
        self.integrator2 = 0;


    def propagateStates(self, dt):
        # propagate full nonlinear equations of motion

        # NOTE: This uses "Euler's First-Order Forward Method" which is the
        # numerical integration equivalent of the forward difference formula.
        # We can look into a more sophisticated approach such as Runge-Kutta.
        # Matlab's "ode45" does Runge-Kutta with 4th and 5th order pair.

        # use states, inputs, and constants to find the derivatives of states
        # first find the drag coefficient based on the current angle theta
        # then use the change in height to find the current air density
        # then find the overall drag and the subsequent acceleration

        # simulate the transient behavior of the air brake
        #self.airBrake(dt)

        CD      = self.CD_b + self.CD_s * self.th;

        h_diff  = self.h - self.h_0
        ratio   = (self.T_0 - self.alpha*h_diff)/self.T_0
        rho     = self.rho_0 * ratio**(self.n-1)

        drag    = 0.5*rho*(self.hd**2)*CD*self.A_ref
        hdd     = -self.g - (drag/self.m)*np.sign(self.hd)

        thdd    = self.u - 0.5*rho*(self.hd**2)*(self.CD_s * self.th)*self.A_ref

        # update the states
        self.h  = self.h  + dt*self.hd
        self.hd = self.hd + dt*hdd

        self.th = self.th + dt*self.thd
        self.thd= self.thd+ dt*thdd

        # populate the state truth vector
        self.x_tru[0,0] = self.h   # altitude truth
        self.x_tru[1,0] = self.hd  # altitude change truth
        self.x_tru[2,0] = self.th  # theta truth
        self.x_tru[3,0] = self.thd # theta change truth
        self.x_tru[4,0] = self.CD_b*self.A_ref
        self.x_tru[5,0] = self.CD_s
        self.x_tru[6,0] = self.g

        # estimate states
        self.estimateStates(dt)

        # update the history vectors
        self.h_all[self.i]  = self.h
        self.hd_all[self.i] = self.hd
        self.t_all[self.i]  = self.t_all[self.i-1]+dt
        self.th_all[self.i] = self.th
        self.x_tru_all[:,self.i] = self.x_tru.ravel()
        self.x_hat_all[:,self.i] = self.x_hat.ravel()

        # interpolate the current hd_ref and compare to reference trajectory
        self.hd_cmd         = np.interp(self.h, self.h_ref, self.hd_ref)
        self.e_hd[self.i]   = self.hd - self.hd_cmd



        # increment index
        self.i = self.i + 1

    def airBrake(self, dt):
        # this simulates the constant-change model of the air brake
        # using the commanded and current theta (th_cmd & th)

        # determine the magnitude of the change in angle
        d_theta_mag = min(self.th_r*dt, abs(self.th_cmd - self.th))

        # determine the direction of the change in angle
        d_theta_sign = np.sign(self.th_cmd - self.th)

        # propagate the air brake angle and update the history vector
        self.th = self.th + d_theta_mag*d_theta_sign
        self.th_all[self.i] = self.th

    def refTrajectory(self, CD_ref):
        # use pessimistic parameters to define a reference trajectory (hd vs h)
        # this is done starting from sea level in order to be more complete
        # vary hd_0_ref until the open loop trajectory arrives at desired h

        # in the future, if we have a good on-line drag estimation, this
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
            h_ref   = np.array([h])
            hd_ref  = np.array([hd])

            # run the trajectory with the current parameters
            while hd > 0:
                # find air density at this altitude
                h_diff  = h - self.h_0
                ratio   = (self.T_0 - self.alpha*h_diff)/self.T_0
                rho     = self.rho_0 * ratio**(self.n-1)

                # get drag and deceleration
                drag    = 0.5*rho*(hd**2)*CD_ref*self.A_ref
                hdd     = -self.g - (drag/self.m)*np.sign(hd)

                # update the states with crude foward difference
                # could be improved with numerical method
                h  = h  + dt*hd
                hd = hd + dt*hdd

                # update the history vectors
                # the if check is to avoid writing early, saving time
                if disparity < thresh*10:
                    h_ref  = np.append(h_ref,  [h] )
                    hd_ref = np.append(hd_ref, [hd])

            # update the inital hd based on arrival trajectory
            disparity = self.h_f - h
            hd_0_ref = hd_0_ref + 1.5*(disparity)

        # save the final reference trajectory as a class member
        self.h_ref  = h_ref
        self.hd_ref = hd_ref

        # generate a starting hd_cmd and the first error history
        self.hd_cmd     = np.interp(self.h, self.h_ref, self.hd_ref)
        self.e_hd[0]    = self.hd - self.hd_cmd
        self.error_d1   = self.hd - self.hd_cmd
        self.error2_d1  = self.th - self.th_cmd

    def saturateControl(self):
        # prevent the control input from going beyond physical constraints
        self.th_cmd = max(self.th_cmd, 0)
        self.th_cmd = min(self.th_cmd, self.th_max)

    def setControl(self, dt):
        # this is the basic controller for now

        # add logic to check if we've passed the "safe deploy velocity"

        # sim is currently using a guessed model for the relationship between
        # theta and DC

        error = self.hd - self.hd_cmd
        error_dot = (error - self.error_d1)/dt
        self.integrator = self.integrator + (self.error_d1 + error)*dt/2

        proportional = error        * 0.5
        integral = self.integrator  * 0.2
        derivative = error_dot      * 0   # no derivative right now

        self.th_cmd = proportional + derivative + integral
        self.saturateControl()

        # add anti-windup logic

        # remember this iteration's error
        self.error_d1 = error

        # remember the control input for plotting
        self.th_cmd_all[self.i] = self.th_cmd


    def setControl2(self, dt):
        # this is the basic controller for dynamic theta

        error = self.th_cmd - self.th
        error_dot = (error - self.error2_d1)/dt
        self.integrator2 = self.integrator2 + (self.error2_d1 + error)*dt/2

        proportional = error        * 10
        integral = self.integrator  * 0
        derivative = error_dot      * 4

        self.u = proportional + derivative + integral
        #self.saturateControl()

        # add anti-windup logic

        # remember this iteration's error
        self.error2_d1 = error

        # remember the control input for plotting
        #self.th_cmd_all[self.i] = self.th_cmd


    def estimateStates(self, dt):
        # this is the estimation
        # simulate noise on the sensors and

        # extract states for ease of coding
        x1 = self.x_hat[0,0]
        x2 = self.x_hat[1,0]
        x3 = self.x_hat[2,0]
        x4 = self.x_hat[3,0]
        x5 = self.x_hat[4,0]
        x6 = self.x_hat[5,0]
        x7 = self.x_hat[6,0]

        # populate A matrix
        frac = (self.T_0 - self.alpha*(x1-self.h_0))/self.T_0
        term = self.rho_0 * frac**(self.n-1)

        front = self.rho_0 * self.alpha * x2**2 / (2 * self.T_0)
        self.A[1,0] = front * (x5 + x6*x3) * (self.n-1) * frac**(self.n-2)
        self.A[1,1] = term * x2 * (x5 + x6*x3) * -1
        self.A[1,2] = term * x2**2 * x6 * -0.5
        self.A[1,4] = term * x2**2 * -0.5
        self.A[1,5] = term * x2**2 * x3 * -0.5

        self.A[3,0] = front * x6 * x3 * (self.n-1) * frac**(self.n-2)
        self.A[3,1] = term * x2 * x6 * x3 * -1
        self.A[3,2] = term * x2**2 * x6 * -0.5
        self.A[3,5] = term * x2**2 * x3 * -0.5

        # generate state transition matrix
        F = linalg.expm(self.A*dt)


        #x_hat_dot = A*self.x_hat
        #self.x_hat = self.x_hat + x_hat_dot*dt

        #P =
        #

        # get sensor data with truth + noise
        y_meas = np.zeros((3,1))

        h_var = 10 # variance on altitude measurement MAKE THESE CLASS MEMBERS FOR THE R MATRIX
        th_var = 1 # variance on theta measurement

        y_meas[0,0] = self.h   + np.random.normal(0, np.sqrt(h_var))
        y_meas[1,0] = self.th  + np.random.normal(0, np.sqrt(th_var))
        y_meas[2,0] = self.g
        #print(y_meas)

        # innovation term
        #C =
        #y_hat = C*x_hat


        # generate kalman gain

        # perform update


        #print('here')
        #print(self.x_hat)
