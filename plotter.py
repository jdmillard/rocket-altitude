from pyqtgraph.Qt import QtGui, QtCore
import numpy as np
import pyqtgraph as pg
import time

class livePlotter:
    """
    Class for plotting methods.
    """
    def __init__(self, final_time, plot_real_time):
        # store some inputs
        self.plot_real_time = plot_real_time
        self.tf = final_time

        ''' setup real time plot using pyqtgraph '''
        if self.plot_real_time:
            # set an application handle
            self.app = QtGui.QApplication([])

            # create the widget ("Graphics Window" allows stacked plots)
            self.win = pg.GraphicsWindow(title="Live Plotting")
            self.win.resize(1500,500)    # set window size
            self.win.move(50,50)         # set window monitor position
            self.win.setWindowTitle('Altitude Controller')

            # enable antialiasing for prettier plots
            pg.setConfigOptions(antialias=True)

            # set some pen types
            pen_green = pg.mkPen(color=(50, 255, 50, 255), width=4)
            pen_blue = pg.mkPen(color=(50, 50, 255, 255), width=2, symbol='t')
            pen_blue2 = pg.mkPen(color=(50, 50, 255, 255), width=1)

            # FIRST SUBPLOT OBJECT
            self.p1 = self.win.addPlot(title="Altitude vs. Time")
            self.p1.setXRange(0,final_time,padding=0)
            self.p1.setYRange(1100,3048,padding=0)
            self.p1.setLabel('left', "Altitude (m)")
            self.p1.setLabel('bottom', "Time (s)") # , units='s'
            self.p1.showGrid(x=True, y=True)
            self.meas1 = self.p1.plot(pen=pen_blue, name='Traj')
            #self.meas1_fade =  self.p1.plot(pen=None, symbol='o', symbolPen=None, symbolSize=6, symbolBrush=(255, 255, 51, 30))

            # SECOND SUBPLOT OBJECT
            self.p2 = self.win.addPlot(title="h_dot vs. h")
            self.p2.setXRange(1100,3048,padding=0)
            self.p2.setYRange(0,320,padding=0)
            self.p2.setLabel('left', "h_dot (m/s)")
            self.p2.setLabel('bottom', "h (m)")
            self.p2.showGrid(x=True, y=True)
            self.meas2 = self.p2.plot(pen=pen_blue, name='Traj2')
            #self.meas2 = self.p2.plot(pen=None, symbol='o', symbolPen=None, symbolSize=3, symbolBrush=(255, 100, 0, 150))
            #self.meas2_fade =  self.p2.plot(pen=None, symbol='o', symbolPen=None, symbolSize=6, symbolBrush=(255, 100, 0, 30))

            '''
            # THIRD SUBPLOT OBJECT
            self.p3 = self.win.addPlot(title="Tracks")
            self.p3.setXRange(-11,11,padding=0)
            self.p3.setYRange(-11,11,padding=0)
            #self.p3.setLabel('left', "Y Axis", units='A')
            self.p3.setLabel('bottom', "East")
            self.p3.showGrid(x=True, y=True)
            # add a legend to subplot, note: you could add handle, then use
            # .addItem, but using name='' in plot objects adds them for you
            self.p3.addLegend()
            # plot the truth trajectory (no data specified yet)
            # the handle will be used with .setData(,) to update live
            self.truth = self.p3.plot(pen=pen_green, name='Truth')
            self.estim = self.p3.plot(pen=pen_blue, name='Estimate')
            self.ellip = self.p3.plot(pen=pen_blue2)
            '''


            # show the plot by calling an update
            # it is needed twice (to force display on first iteration) - not sure why
            # either method below works, but the app handle method is better practice
            self.app.processEvents() #pg.QtGui.QApplication.processEvents()
            self.app.processEvents() #pg.QtGui.QApplication.processEvents()

            # start timer
            self.time0 = time.time()


    # method for updating data
    def updateItems(self, rocket, sim_time, current_time):
        if self.plot_real_time:
            # plot no faster than actual time
            actual_time = current_time - self.time0
            if actual_time < sim_time:
                # pause to wait for actual time to catch up
                time.sleep(sim_time-actual_time)

            # get time and h for the rocket
            x = rocket.t_all[0:rocket.i]
            y = rocket.h_all[0:rocket.i]
            self.meas1.setData(x,y)

            # get h and h_dot for the rocket
            x = rocket.h_all[0:rocket.i]
            y = rocket.hd_all[0:rocket.i]
            self.meas2.setData(x,y)

            #x = pda.sensor_list[0].meas_list[:,0]
            #y = pda.sensor_list[0].meas_list[:,1]
            #self.meas1_fade.setData(x,y)

            '''
            # get most recent measurement from sensor 2
            x = pda.sensor_list[1].meas[:,0]
            y = pda.sensor_list[1].meas[:,1]
            self.meas2.setData(x,y)
            x = pda.sensor_list[1].meas_list[:,0]
            y = pda.sensor_list[1].meas_list[:,1]
            self.meas2_fade.setData(x,y)

            # pull out east and north elements
            x = agents.state_list[:,1]
            y = agents.state_list[:,0]
            self.truth.setData(x,y)
            x = pda.x_hat_list[:,0]
            y = pda.x_hat_list[:,1]
            self.estim.setData(x, y)

            # get first two states and corresponding covariance
            in1 = pda.x_hat[0:2,:]
            in2 = pda.P[0:2,0:2]

            increment = 0.1
            phi = np.arange(0, 2*np.pi+increment, increment)
            aa = np.atleast_2d(np.cos(phi))
            bb = np.atleast_2d(np.sin(phi))
            pts = np.concatenate((aa,bb))

            U = np.linalg.cholesky(in2)
            sigma = 3

            pts = sigma*np.matmul(np.transpose(U),pts) + np.tile(in1, len(phi))
            self.ellip.setData(pts[0,:],pts[1,:])
            '''

            # update the plotted data
            self.app.processEvents() #pg.QtGui.QApplication.processEvents()

            # hold plot on last iteration
            print('---')
            print(sim_time)
            print(self.tf)
            if int(sim_time) == self.tf:
                print("desired simulation time: ", self.tf, ", time taken: ", current_time - self.time0)
                self.app.exec_() # hold final plot

    # method for generating 2d ellipse for a given covariance
    def generateEllipse(self, P):
        # fill in ellipse generation here
        return 3
