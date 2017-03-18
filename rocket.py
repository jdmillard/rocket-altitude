import numpy as np

# this early in the repo history, this is all artifacts of what I pasted
# from another project

class PDAfilter:
    """
    Class for PDA filter. A given PDA filter is assigned to 1 target.
    With imf, there are multiple sensors, so a list of noise sigmas is sent.
    """
    def __init__(self, dt, s_list, agent):
        # the length of the s list determines the number of sensors
        # the sensor list contains the objects for each sensor
        # agent object is passed in for initial sensor measurement
        # initialize them here:
        self.sensor_list = []
        for params in s_list:
            self.sensor_list.append(sensor(params, agent))

        # the structure of the filter accommodates variable dt
        # t_last is for finding elapsed time between propagations
        # t_last_u is for finding the elapsed time between measurements
        self.t_last = 0
        self.t_last_u = 0

        # x_hat is current estimate, x_hat_u is the last update estimate
        # x_hat can be propagated forward and x_hat_u remembered
        self.x_hat = np.zeros((4,1))
        self.x_hat[0] = agent.state[1]
        self.x_hat[1] = agent.state[0]
        self.x_hat[2] = np.pi           # MANUALLY ENTERED INITIAL VELOCITY
        self.x_hat_u = self.x_hat

        # save state as a row at each time step
        self.x_hat_list = np.transpose(self.x_hat)

        # initialize P based on confidence of initial state
        # these are initialized at truth, hence low values: x, y, x_dot, y_dot
        self.P = np.zeros((4,4))
        self.P[0,0] = 0.01
        self.P[1,1] = 0.01
        self.P[2,2] = 0.01
        self.P[3,3] = 0.01

        # process noise; this can be tuned
        self.Q = np.zeros((4,4))
        self.Q[0,0] = 0.1
        self.Q[1,1] = 0.1
        self.Q[2,2] = 1.5
        self.Q[3,3] = 1.5

        self.P_information = np.linalg.inv(self.P)
        self.x_information = np.matmul(self.P_information, self.x_hat)


    def getUpdate(self, agent, t):
        # have each sensor perform measurement
        # target object is passed in order to use source truth
        for sensor in self.sensor_list:
            sensor.getMeasurement(agent, t)


    def propagateStates(self, t):
        # if t==0, the state and estimates are already caught up to time
        # and therefore have no need of propagation
        if t==0:
            self.t_last = t
            return

        # generate time since last propagation or measurement update
        t_elapsed = t - self.t_last
        self.t_last = t

        # generate the state transition matrix based on elapsed time
        F = np.zeros((4,4))
        F[0,0] = 1
        F[0,2] = t_elapsed
        F[1,1] = 1
        F[1,3] = t_elapsed
        F[2,2] = 1
        F[3,3] = 1

        # update current state estimate and estimate covariance
        F_t = np.transpose(F)
        self.x_hat = np.matmul(F, self.x_hat)
        self.P = np.matmul(np.matmul(F, self.P), F_t) + self.Q * t_elapsed

        # NOTE: in an EKF setting, x_hat could be added to x_hat_list
        # currently, it's all done at the same time so it would only be clutter


    def executeFilter(self, t):
        # if t==0, the state and estimates are already caught up to time
        # and therefore have no need of propagation
        if t == 0:
            self.t_last_u = t
            return

        # generate time since last measurement update
        # this is used for crude velocity measurements between points
        t_elapsed_u = t - self.t_last_u
        self.t_last_u = t



        # update the post-propagation information matrix and vector
        # then reset the information gain sums before the loop
        self.P_information = np.linalg.inv(self.P)                              # CAN BE MOVED TO PROPAGATION
        self.x_information = np.matmul(self.P_information, self.x_hat)          # CAN BE MOVED TO PROPAGATION
        inf_mat_sum = np.zeros((4,4))
        inf_vec_sum = np.zeros((4,1))
        # execute the PDA filters associated with each sensor, consider the
        # contribution of each one that had measurement updates
        # NOTE: modify structure such that each PDAF can be tuned independently
        for sensor in self.sensor_list:
            x_out, P_out = sensor.filterContribute(t_elapsed_u, self.x_hat, self.P, self.x_hat_u)

            # calculation of "information gain" (increase) for each sensor
            # for a sensor with no measurements, the gain will be zero
            # first convert the x,P from the current sensor into information
            inf_mat_i = np.linalg.inv(P_out)
            inf_vec_i = np.matmul(inf_mat_i, x_out)
            inf_mat_sum = inf_mat_sum + (inf_mat_i - self.P_information)
            inf_vec_sum = inf_vec_sum + (inf_vec_i - self.x_information)


        self.P_information = self.P_information + inf_mat_sum
        self.x_information = self.x_information + inf_vec_sum

        self.P = np.linalg.inv(self.P_information)
        self.x_hat = np.matmul(self.P, self.x_information)

        # CURRENT ISSUE: first few are validated, after it diverges,
        # covariance is not growing - RED FLAG, maybe becaue P_information needs to be updated post-propagation ???


        # hacky way of only using the last PDA filter
        #self.x_hat = x_out
        #self.P     = P_out

        # update the remembered estimate associated with the most recent update
        self.x_hat_u = self.x_hat

        # update history of estimates
        self.x_hat_list = np.concatenate((self.x_hat_list, np.transpose(self.x_hat)), axis=0)




class sensor:
    """
    Class for sensor
    """
    def __init__(self, params, agent):
        self.params = params
        self.noise = params['noise']      # standard deviation
        self.R = np.zeros((4,4))
        self.R[0,0] = self.noise
        self.R[1,1] = self.noise
        self.R[2,2] = self.noise
        self.R[3,3] = self.noise
        self.n_meas = params['n_meas']

        # generate an initial measurement
        self.meas = np.zeros(2)
        self.meas[0] = agent.state[1] + np.random.normal(0, self.noise, 1)
        self.meas[1] = agent.state[0] + np.random.normal(0, self.noise, 1)

        # measurement history
        self.meas_list = np.atleast_2d(self.meas)

        # a gamma value is chosen as a gate threshold for the validation region
        # based on the Mahalanobis distance of point measurements.
        # (this is sigma value squared; gamma = 9 for 3 sigma validation)
        # it is okay to validate 3 sigma because association likelihood values
        # are weighted by multivariate normal distribution
        self.H = np.eye(4)  # simple sensor measurement model
        self.gamma = 9      # 3 sigma validation region
        self.PD    = 0.9    # detection probability
        self.PG    = 0.9    # gate probability
        self.lambd = 0.01   # poisson spatial density of interference points


    def getMeasurement(self, agent, t):
        # simulate the measurements for each sensor
        if t<self.params['dead'] or t>self.params['alive']:
            self.meas = np.zeros((self.n_meas,2))
            self.meas[0:self.n_meas,0] = agent.state[1] + np.random.normal(0, self.noise, self.n_meas)
            self.meas[0:self.n_meas,1] = agent.state[0] + np.random.normal(0, self.noise, self.n_meas)

            # update measurement history
            self.meas_list = np.concatenate((self.meas_list, self.meas), axis=0)
        else:
            self.meas = np.ones((1,2))*100


    def filterContribute(self, dt_u, x_hat, P, x_hat_u):
        # this method performs a PDA measurement update for the given sensor.
        # the input is the result of the previous measurements' IMF after
        # recent propagation. the output from all sensors goes to IMF

        # calculate expected measurement and innovation covariance
        z_hat = np.matmul(self.H, x_hat) # standard EKF adds sample N(0,R_k)
        H_t = np.transpose(self.H)
        S = np.matmul(np.matmul(self.H, P), H_t) + self.R
        S_inv = np.linalg.inv(S) # NOTE: look into chol instead of inverse

        # for each measurement (each row), check against validation region
        # CURRENTLY, ONLY 1 MEASUREMENT PER SENSOR, so force atleast_2d
        # if interference is added, z_batch will come 2d
        # L is likelihood vecor (one for each measurement)
        z_batch = self.meas
        L = np.zeros((z_batch.shape[0], 1))
        for ind_m in range(z_batch.shape[0]):
            z_cur = np.zeros((1,4))             # reset z_cur
            z_cur[0,0:2] = z_batch[ind_m, :]    # get x, y positions
            z_cur = np.transpose(z_cur)         # set upright

            # create velocity "measurement" using last update values
            z_cur[2,0] = (z_cur[0,0]-x_hat_u[0,0])/dt_u # x_dot
            z_cur[3,0] = (z_cur[1,0]-x_hat_u[1,0])/dt_u # y_dot

            # individual innovation vector
            nu_i = z_cur - z_hat
            nu_i_t = np.transpose(nu_i)

            # use nu_i and S to check elliptical validation region
            # nu_i determines distance from "origin"
            d_mah = np.matmul(np.matmul(nu_i_t, S_inv), nu_i)[0][0]

            if d_mah <= self.gamma:
                # current measurement is validated
                # get association likelihood (Li) using multivariate gaussian

                # get non-normalized surface curve value
                # then normalize
                # consider detection probability and poisson interference
                curve = np.exp(-0.5*d_mah)
                dim = 4 # state dimensionality, MANUALLY ENTERED
                Li = curve/np.sqrt(np.linalg.det(S)*np.power(2*np.pi,dim/2))
                Li = Li * self.PD / self.lambd

                # store association likelihood values for each measurement
                # unvalidated measurements will default to likelihood of 0
                L[ind_m,0]=Li

        # now L contains the likelihoods of all measurements
        # convert it to association probability
        # this is more complicated than normalizing the vector because we
        # need to account for the probability of no received measurement
        # if PD=PG=1, the following breaks down to standard normalization
        beta = L / (1 - self.PD*self.PG + np.sum(L))

        # probability that none are associated, not part of vector
        beta_0 = (1 - self.PD*self.PG) / (1 - self.PD*self.PG + np.sum(L))

        # find nu innovation vector (weighted sum of all individual nu_i)
        # cannot be merged with above loop because we first need beta vector
        # also, get 'spread' term for the summation needed below (see 3.5.2-16)
        nu = np.zeros((4,1))
        spread = np.zeros((4,1))
        for ind_m in range(z_batch.shape[0]):
            z_cur = np.zeros((1,4))             # reset z_cur
            z_cur[0,0:2] = z_batch[ind_m, :]    # get x, y positions
            z_cur = np.transpose(z_cur)         # set upright

            # create velocity "measurement" using last update values
            z_cur[2,0] = (z_cur[0,0]-x_hat_u[0,0])/dt_u # x_dot
            z_cur[3,0] = (z_cur[1,0]-x_hat_u[1,0])/dt_u # y_dot

            # individual innovation vector
            nu_i = z_cur - z_hat
            nu_i_t = np.transpose(nu_i)

            nu = nu + (nu_i*beta[ind_m,0])
            spread = spread + (beta[ind_m, 0] * np.matmul(nu_i, nu_i_t))

        # calculate the gain
        W = np.matmul(np.matmul(P, H_t), S_inv)

        # update state estimate
        x_hat = x_hat + np.matmul(W, nu)

        # calculate the covariance associated with the correct measurement
        W_t = np.transpose(W)
        Pc = P - np.matmul(np.matmul(W, S), W_t)

        # calculate the spread of the innovations term (see 3.5.2-16)
        nu_t = np.transpose(nu)
        Ptilde = np.matmul(np.matmul(W, (spread - np.matmul(nu, nu_t))), W_t)

        # calculate the covariance of the updated state
        P = beta_0*P + (1-beta_0)*Pc + Ptilde

        # update complete
        return x_hat, P
