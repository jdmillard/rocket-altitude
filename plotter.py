from pyqtgraph.Qt import QtGui, QtCore
import numpy as np
import pyqtgraph as pg
import time

class livePlotter:
    """
    Class for plotting methods.
    """
    def __init__(self, rocket, final_time, plot_real_time):
        # store some inputs
        self.plot_real_time = plot_real_time
        self.tf = final_time

        ''' setup real time plot using pyqtgraph '''
        self.app = QtGui.QApplication([])

        # create the widget ("Graphics Window" allows stacked plots)
        self.win0 = pg.GraphicsWindow(title="Live Plotting")
        self.win0.resize(1500,1000)    # set window size
        self.win0.move(50,50)          # set window monitor position
        self.win0.setWindowTitle('Altitude Controller Truth')

        # create another widget ("Graphics Window" allows stacked plots)
        self.win1 = pg.GraphicsWindow(title="Live Plotting")
        self.win1.resize(1500,600)    # set window size
        self.win1.move(200,200)          # set window monitor position
        self.win1.setWindowTitle('Estimation')

        # enable antialiasing for prettier plots
        pg.setConfigOptions(antialias=True)

        # set some pen types
        pen_green  = pg.mkPen(color=(50, 255, 50, 255), width=2)
        pen_green2 = pg.mkPen(color=(50, 255, 50, 255), width=1)
        pen_blue   = pg.mkPen(color=(50, 50, 255, 255), width=2, symbol='t')
        pen_blue2  = pg.mkPen(color=(50, 50, 255, 255), width=1)

        # FIRST SUBPLOT OBJECT
        self.p1 = self.win0.addPlot(title="Altitude vs. Time")
        self.p1.setXRange(0,final_time,padding=0)
        self.p1.setYRange(rocket.h*0.9,rocket.h_f*1.1,padding=0)
        self.p1.setLabel('left', "Altitude (m)")
        self.p1.setLabel('bottom', "Time (s)") # , units='s'
        self.p1.showGrid(x=True, y=True)
        self.meas1 = self.p1.plot(pen=pen_blue, name='Curve 1')

        # SECOND SUBPLOT OBJECT
        self.p2 = self.win0.addPlot(title="Velocity vs. Time")
        self.p2.setXRange(0,final_time,padding=0)
        self.p2.setYRange(0,rocket.hd_0*1.1,padding=0)
        self.p2.setLabel('left', "h_dot (m/s)")
        self.p2.setLabel('bottom', "Time (s)")
        self.p2.showGrid(x=True, y=True)
        self.meas2 = self.p2.plot(pen=pen_blue, name='Curve 2')

        # THIRD SUBPLOT OBJECT
        self.p3 = self.win0.addPlot(title="h_dot vs. h")
        self.p3.setXRange(rocket.h*0.9,rocket.h_f*1.1,padding=0)
        self.p3.setYRange(0,rocket.hd_0*1.1,padding=0)
        self.p3.setLabel('left', "h_dot (m/s)")
        self.p3.setLabel('bottom', "h (m)")
        self.p3.showGrid(x=True, y=True)
        self.p3.addLegend(offset=[-10,10])
        self.meas3 = self.p3.plot(pen=pen_blue, name='Simulated Trajectory')
        self.t_ref = self.p3.plot(pen=pen_green2, name='Reference Trajectory')
        self.t_ref.setData(rocket.h_ref, rocket.hd_ref)

        self.win0.nextRow()

        # FOURTH SUBPLOT OBJECT
        self.p4 = self.win0.addPlot(title="Theta Control Input")
        self.p4.setXRange(0,final_time,padding=0)
        self.p4.setYRange(0,rocket.th_max*1.1,padding=0)
        self.p4.setLabel('left', "theta (deg)")
        self.p4.setLabel('bottom', "time (s)")
        self.p4.showGrid(x=True, y=True)
        self.p4.addLegend(offset=[-10,10])
        self.meas4 = self.p4.plot(pen=pen_blue, name='Current Theta')
        self.meas4a = self.p4.plot(pen=pen_green2, name='Desired Theta')

        # FIFTH SUBPLOT OBJECT
        self.p5 = self.win0.addPlot(title="Error vs. Time")
        #self.p5.setLogMode(False,True)
        self.p5.setXRange(0,final_time,padding=0)
        #self.p5.setYRange( , ,padding=0)
        self.p5.setLabel('left', "Velocity Error (m/s)")
        self.p5.setLabel('bottom', "Time (s)")
        self.p5.showGrid(x=True, y=True)
        self.meas5 = self.p5.plot(pen=pen_green, name='Curve 6')

        # SIXTH SUBPLOT OBJECT
        self.p6 = self.win0.addPlot(title="Error vs. Height")
        self.p6.setXRange(rocket.h*0.9,rocket.h_f*1.1,padding=0)
        #self.p6.setYRange(rocket.h*0.9,rocket.h_f*1.1,padding=0)
        self.p6.setLabel('left', "Velocity Error (m/s)")
        self.p6.setLabel('bottom', "h (m)")
        self.p6.showGrid(x=True, y=True)
        self.meas6 = self.p6.plot(pen=pen_green, name='Curve 6')

        # SECOND WINDOW

        # FIRST SUBPLOT OBJECT
        self.q1 = self.win1.addPlot(title="h")
        #self.q1.setXRange(0,final_time,padding=0)
        #self.q1.setYRange(rocket.h*0.9,rocket.h_f*1.1,padding=0)
        self.q1.setLabel('left', "Altitude (m)")
        self.q1.setLabel('bottom', "Time (s)") # , units='s'
        self.q1.showGrid(x=True, y=True)
        self.state01_t = self.q1.plot(pen=pen_green, name='Truth')
        self.state01_e = self.q1.plot(pen=pen_blue, name='Estimate')

        # SECOND SUBPLOT OBJECT
        #self.win1.nextRow()
        self.q2 = self.win1.addPlot(title="h_dot")
        #self.q2.setXRange(0,final_time,padding=0)
        #self.q2.setYRange(rocket.h*0.9,rocket.h_f*1.1,padding=0)
        self.q2.setLabel('left', "Altitude Change (m/s)")
        self.q2.setLabel('bottom', "Time (s)") # , units='s'
        self.q2.showGrid(x=True, y=True)
        self.state02_t = self.q2.plot(pen=pen_green, name='Truth')
        self.state02_e = self.q2.plot(pen=pen_blue, name='Estimate')

        # THIRD SUBPLOT OBJECT
        self.win1.nextRow()
        self.q3 = self.win1.addPlot(title="theta")
        #self.q3.setXRange(0,final_time,padding=0)
        #self.q3.setYRange(rocket.h*0.9,rocket.h_f*1.1,padding=0)
        self.q3.setLabel('left', "Theta (deg)")
        self.q3.setLabel('bottom', "Time (s)") # , units='s'
        self.q3.showGrid(x=True, y=True)
        self.state03_t = self.q3.plot(pen=pen_green, name='Truth')
        self.state03_e = self.q3.plot(pen=pen_blue, name='Estimate')

        # FOURTH SUBPLOT OBJECT
        #self.win1.nextRow()
        self.q4 = self.win1.addPlot(title="theta dot")
        #self.q4.setXRange(0,final_time,padding=0)
        #self.q4.setYRange(rocket.h*0.9,rocket.h_f*1.1,padding=0)
        self.q4.setLabel('left', "Theta Dot (deg/s)")
        self.q4.setLabel('bottom', "Time (s)") # , units='s'
        self.q4.showGrid(x=True, y=True)
        self.state04_t = self.q4.plot(pen=pen_green, name='Truth')
        self.state04_e = self.q4.plot(pen=pen_blue, name='Estimate')

        # FIFTH SUBPLOT OBJECT
        self.win1.nextRow()
        self.q5 = self.win1.addPlot(title="CD_0*A")
        #self.q5.setXRange(0,final_time,padding=0)
        #self.q5.setYRange(rocket.h*0.9,rocket.h_f*1.1,padding=0)
        self.q5.setLabel('left', "Base Drag")
        self.q5.setLabel('bottom', "Time (s)") # , units='s'
        self.q5.showGrid(x=True, y=True)
        self.state05_t = self.q5.plot(pen=pen_green, name='Truth')
        self.state05_e = self.q5.plot(pen=pen_blue, name='Estimate')

        # SIXTH SUBPLOT OBJECT
        #self.win1.nextRow()
        self.q6 = self.win1.addPlot(title="CD_0*A change")
        #self.q6.setXRange(0,final_time,padding=0)
        #self.q6.setYRange(rocket.h*0.9,rocket.h_f*1.1,padding=0)
        self.q6.setLabel('left', "Drag Rate")
        self.q6.setLabel('bottom', "Time (s)") # , units='s'
        self.q6.showGrid(x=True, y=True)
        self.state06_t = self.q6.plot(pen=pen_green, name='Truth')
        self.state06_e = self.q6.plot(pen=pen_blue, name='Estimate')

        # show the plot by calling an update
        # it is needed twice (to force display on first iteration) - not sure why
        # either method below works, but the app handle method is better practice
        self.app.processEvents() #pg.QtGui.QApplication.processEvents()
        self.app.processEvents() #pg.QtGui.QApplication.processEvents()

        # start timer
        self.time0 = time.time()


    # method for updating data
    def updateItems(self, rocket, sim_time, current_time):

        # override the waiting constraint
        if self.plot_real_time:
            actual_time = current_time - self.time0
        else:
            actual_time = sim_time


        if self.plot_real_time or rocket.hd <= 0 or sim_time==self.tf:
            # plot no faster than actual time
            # NOTE: simulation can get slower than real time
            if actual_time < sim_time:
                # pause to wait for actual time to catch up
                time.sleep(sim_time-actual_time)


            # get time and h for the rocket
            x = rocket.t_all[0:rocket.i]
            y = rocket.h_all[0:rocket.i]
            self.meas1.setData(x,y)

            # get time and h_dot for the rocket
            #x = rocket.t_all[0:rocket.i] # x is already this
            y = rocket.hd_all[0:rocket.i]
            self.meas2.setData(x,y)

            # get h and h_dot for the rocket
            x = rocket.h_all[0:rocket.i]
            #y = rocket.hd_all[0:rocket.i] # y is already this
            self.meas3.setData(x,y)

            # get time and theta for the air brake
            x = rocket.t_all[0:rocket.i]
            y = rocket.th_all[0:rocket.i]
            self.meas4.setData(x,y)

            # get time and theta_cmd for the air brake
            #x = rocket.t_all[0:rocket.i]
            y = rocket.th_cmd_all[0:rocket.i]
            self.meas4a.setData(x,y)

            # get time and e_hd for the rocket
            #x = rocket.t_all[0:rocket.i]
            y = rocket.e_hd[0:rocket.i]
            self.meas5.setData(x,y)

            # get h and e_hd for the rocket
            x = rocket.h_all[0:rocket.i]
            #y = rocket.e_hd[0:rocket.i]
            self.meas6.setData(x,y)

            # get h truth and estimate for the rocket
            x = rocket.t_all[0:rocket.i]
            y = rocket.x_tru_all[0,0:rocket.i]
            self.state01_t.setData(x,y)
            y = rocket.x_hat_all[0,0:rocket.i]
            self.state01_e.setData(x,y)

            # get h_dot truth and estimate for the rocket
            #x = rocket.t_all[0:rocket.i]
            y = rocket.x_tru_all[1,0:rocket.i]
            self.state02_t.setData(x,y)
            y = rocket.x_hat_all[1,0:rocket.i]
            self.state02_e.setData(x,y)

            # get theta truth and estimate for the rocket
            #x = rocket.t_all[0:rocket.i]
            y = rocket.x_tru_all[2,0:rocket.i]
            self.state03_t.setData(x,y)
            y = rocket.x_hat_all[2,0:rocket.i]
            self.state03_e.setData(x,y)

            # get theta dot truth and estimate for the rocket
            #x = rocket.t_all[0:rocket.i]
            y = rocket.x_tru_all[3,0:rocket.i]
            self.state04_t.setData(x,y)
            y = rocket.x_hat_all[3,0:rocket.i]
            self.state04_e.setData(x,y)

            # get theta dot truth and estimate for the rocket
            #x = rocket.t_all[0:rocket.i]
            y = rocket.x_tru_all[4,0:rocket.i]
            self.state05_t.setData(x,y)
            y = rocket.x_hat_all[4,0:rocket.i]
            self.state05_e.setData(x,y)

            # get theta dot truth and estimate for the rocket
            #x = rocket.t_all[0:rocket.i]
            y = rocket.x_tru_all[5,0:rocket.i]
            self.state06_t.setData(x,y)
            y = rocket.x_hat_all[5,0:rocket.i]
            self.state06_e.setData(x,y)


            # update the plotted data
            self.app.processEvents() #pg.QtGui.QApplication.processEvents()

            # hold plot when rocket reaches maximum height
            if rocket.hd <= 0 or sim_time==self.tf:
                print("simulation finished")
                print("rocket altitude:", rocket.h, "m")
                print("simulation time:", sim_time, "s")
                print(rocket.cum_error)
                #print("real time: ", current_time - self.time0, " s")
                while 1:
                    self.app.processEvents() #pg.QtGui.QApplication.processEvents()
                    self.app.exec_() # hold final plot
                    #time.sleep(5)

    # method for generating 2d ellipse for a given covariance
    def generateEllipse(self, P):
        # fill in ellipse generation here
        return 3
