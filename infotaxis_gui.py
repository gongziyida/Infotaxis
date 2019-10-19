"""
Author: Ziyi Gong
"""
from pyqtgraph.Qt import QtCore, QtGui
from threading import Thread
import pyqtgraph as pg
import numpy as np
import random
import sys
import time
import os

import infotaxis as itx

np.set_printoptions(precision=2)


class InfotaxisRealTime(QtGui.QWidget):
    def __init__(self, infotaxis, freq=0.01):
        """
        Parameters
        ----------
        infotaxis: infotaxis.infotaxis
            An infotaxis object, the "searcher"
        freq: float, optional
            The refreshing frequency for the animation
        """
        super().__init__()

        # the simulation module
        self.infotaxis = infotaxis

        # Desired Frequency (Hz) = 1 / self.FREQUENCY
        # USE FOR TIME.SLEEP (s)
        self.FREQUENCY = freq

        # Frequency to update plot (ms)
        # USE FOR TIMER.TIMER (ms)
        self.TIMER_FREQUENCY = self.FREQUENCY * 1000

        # Outer grid layout
        self.layout = QtGui.QGridLayout()

        self._make()

    def _make(self):
        # Make the environment widget
        self.env_widget = pg.PlotWidget()
        self.env_widget.hideAxis('left')
        self.env_widget.hideAxis('bottom')

        # Data
        src_pos, src_radius = self.infotaxis.src
        self.prob_map = self.infotaxis.prob_map
        self.searcher_pos = np.zeros((1, 2), dtype=int)
        self.searcher_pos[0], searcher_size = self.infotaxis.cur
        self.hist_pos = self.searcher_pos.copy()
        self.hits_record = [False]
        dim = self.infotaxis._dim
        v = self.infotaxis._v
        tau = float(self.infotaxis._tau)
        d = float(self.infotaxis._d)
        r = float(self.infotaxis._r)

        # Estimation of the source position (probability map)
        self.prob_map_heat = pg.ImageItem()
        self.env_widget.addItem(self.prob_map_heat)
        self.prob_map_heat.setImage(self.prob_map)

        # Locate the actual source position
        self.src = pg.ScatterPlotItem(pxMode=False, size=src_radius * 2)
        self.env_widget.addItem(self.src)
        self.src.setData(pos=np.array([src_pos]))

        txt = 'Dimension: {}\tWind: {}\nParticle Lifetime: {:.1f} s\t' + \
              'Diffusivity: {:.1f} m^2/s\tEmission Rate {:.1f} Hz'
        self.winddir = pg.TextItem(text=txt.format(dim, v, tau, d, r))
        self.env_widget.addItem(self.winddir)

        # Trajectory
        self.trajectory = pg.PlotCurveItem()
        self.env_widget.addItem(self.trajectory)
        self.trajectory.setData(x=self.hist_pos[:, 0], y=self.hist_pos[:, 1])

        # hits
        self.hits = pg.ScatterPlotItem()
        self.hits.setBrush(50, 0, 0)
        self.env_widget.addItem(self.hits)

        # Searcher
        self.searcher = pg.ScatterPlotItem(pxMode=False, size=searcher_size * 2)
        self.searcher.setBrush(255, 0, 0)
        self.env_widget.addItem(self.searcher)
        self.searcher.setData(pos=self.searcher_pos)

        self.layout.addWidget(self.env_widget)

    def start(self):
        """
        To start the simulation
        """
        # Read in data using a thread
        self.update_thread = Thread(target=self.read_pos, args=())
        self.update_thread.daemon = True
        self.update_thread.start()

        self.update_timer = QtCore.QTimer()
        self.update_timer.timeout.connect(self.plot_updater)
        self.update_timer.start(self.TIMER_FREQUENCY)

    def read_pos(self):
        frequency = self.FREQUENCY
        time.sleep(0.3)  # Make sure the initial step is visible
        while True:
            # sleep
            time.sleep(frequency)

            ret = self.infotaxis.step()

            if ret != -1:
                self.searcher_pos[0], hit = ret

            self.hist_pos = np.append(self.hist_pos, self.searcher_pos, 0)

            self.hits_record.append(hit)

            # update probability map
            self.prob_map = self.infotaxis.prob_map[::-1]

    def plot_updater(self):
        self.trajectory.setData(x=self.hist_pos[:, 0], y=self.hist_pos[:, 1])
        self.hits.setData(pos=self.hist_pos[self.hits_record])
        self.searcher.setData(pos=self.searcher_pos)
        self.prob_map_heat.setImage(self.prob_map)


"""
An example of how to make the animation
"""
if __name__ == '__main__':
    # Create main application window
    app = QtGui.QApplication([])
    app.setStyle(QtGui.QStyleFactory.create("Cleanlooks"))
    mw = QtGui.QMainWindow()
    mw.setWindowTitle('Infotaxis')

    # Create an infotaxis instance
    infotaxis_example = itx.Infotaxis(
        d=2, r=100, a=1, v=0.5, winddir=3 * np.pi / 2, tau=500, dt=0.1, dim=3,
        pos=np.zeros(2), src_pos=np.full(2, 90), src_radius=2, window=100)

    # Make a GUI
    infotaxis_real_time = InfotaxisRealTime(infotaxis_example)

    # Create and set widget layout
    # Main widget container
    cw = QtGui.QWidget()
    cw.setFixedSize(700, 750)
    ml = QtGui.QGridLayout()
    cw.setLayout(ml)
    mw.setCentralWidget(cw)

    # Can use either to add plot to main layout
    ml.addLayout(infotaxis_real_time.layout, 0, 0)
    mw.show()
    infotaxis_real_time.start()

    # Start Qt event loop unless running in interactive mode or using pyside
    if (sys.flags.interactive != 1) or not hasattr(QtCore, 'PYQT_VERSION'):
        QtGui.QApplication.instance().exec_()
