#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#
# SPDX-License-Identifier: GPL-3.0
#
# GNU Radio Python Flow Graph
# Title: Not titled yet
# Author: japjeet
# GNU Radio version: v3.11.0.0git-664-g5884ec6d

from PyQt5 import Qt
from gnuradio import qtgui
from PyQt5 import QtCore
from gnuradio import analog
from gnuradio import audio
from gnuradio import filter
from gnuradio.filter import firdes
from gnuradio import gr
from gnuradio.fft import window
import sys
import signal
from PyQt5 import Qt
from argparse import ArgumentParser
from gnuradio.eng_arg import eng_float, intx
from gnuradio import eng_notation
from gnuradio import uhd
import time

from PyQt5.QtCore import QTimer  # Add this line for QTimer

class manishkitry(gr.top_block, Qt.QWidget):

    def __init__(self):
        gr.top_block.__init__(self, "Not titled yet", catch_exceptions=True)
        Qt.QWidget.__init__(self)
        self.setWindowTitle("Not titled yet")
        qtgui.util.check_set_qss()
        try:
            self.setWindowIcon(Qt.QIcon.fromTheme('gnuradio-grc'))
        except BaseException as exc:
            print(f"Qt GUI: Could not set Icon: {str(exc)}", file=sys.stderr)
        self.top_scroll_layout = Qt.QVBoxLayout()
        self.setLayout(self.top_scroll_layout)
        self.top_scroll = Qt.QScrollArea()
        self.top_scroll.setFrameStyle(Qt.QFrame.NoFrame)
        self.top_scroll_layout.addWidget(self.top_scroll)
        self.top_scroll.setWidgetResizable(True)
        self.top_widget = Qt.QWidget()
        self.top_scroll.setWidget(self.top_widget)
        self.top_layout = Qt.QVBoxLayout(self.top_widget)
        self.top_grid_layout = Qt.QGridLayout()
        self.top_layout.addLayout(self.top_grid_layout)
        self.toggle_frequency = False  # Flag to toggle between two frequencies

        self.timer = QTimer(self)
        self.timer.timeout.connect(self.update_frequency)
        self.timer.start(100)  # Time is in milliseconds (2000 ms = 2 seconds)

        self.settings = Qt.QSettings("GNU Radio", "manishkitry")

        try:
            geometry = self.settings.value("geometry")
            if geometry:
                self.restoreGeometry(geometry)
        except BaseException as exc:
            print(f"Qt GUI: Could not restore geometry: {str(exc)}", file=sys.stderr)

        ##################################################
        # Variables
        ##################################################
        self.volume = volume = 300e-3
        self.tuning = tuning = 98.3e6
        self.sample_rate = sample_rate = 1.92e6
        self.rf_decim = rf_decim = 10
        self.interp = interp = 1
        self.deviation = deviation = 75e3
        self.audio_decim = audio_decim = 4

        ##################################################
        # Blocks
        ##################################################

        self._tuning_range = qtgui.Range(88e6, 108e6, 200e3, 98.3e6, 200)
        self._tuning_win = qtgui.RangeWidget(self._tuning_range, self.set_tuning, "'tuning'", "counter_slider", float, QtCore.Qt.Horizontal)
        self.top_layout.addWidget(self._tuning_win)
        self._volume_range = qtgui.Range(0, 1, 50e-3, 300e-3, 200)
        self._volume_win = qtgui.RangeWidget(self._volume_range, self.set_volume, "'volume'", "counter_slider", float, QtCore.Qt.Horizontal)
        self.top_layout.addWidget(self._volume_win)
        self.uhd_usrp_source_0 = uhd.usrp_source(
            ",".join(("", '')),
            uhd.stream_args(
                cpu_format="fc32",
                args='',
                channels=list(range(0,1)),
            ),
        )
        self.uhd_usrp_source_0.set_samp_rate(sample_rate)
        self.uhd_usrp_source_0.set_time_unknown_pps(uhd.time_spec(0))

        self.uhd_usrp_source_0.set_center_freq(tuning, 0)
        self.uhd_usrp_source_0.set_antenna("TX/RX", 0)
        self.uhd_usrp_source_0.set_bandwidth(500e3, 0)
        self.uhd_usrp_source_0.set_gain(50, 0)
        self.rational_resampler_xxx_0 = filter.rational_resampler_ccc(
                interpolation=interp,
                decimation=10,
                taps=[],
                fractional_bw=0)
        self.audio_sink_0 = audio.sink(48000, '', True)
        self.analog_fm_demod_cf_0 = analog.fm_demod_cf(
        	channel_rate=192e3,
        	audio_decim=audio_decim,
        	deviation=deviation,
        	audio_pass=16000,
        	audio_stop=20000,
        	gain=1.0,
        	tau=(75e-6),
        )


        ##################################################
        # Connections
        ##################################################
        self.connect((self.analog_fm_demod_cf_0, 0), (self.audio_sink_0, 0))
        self.connect((self.rational_resampler_xxx_0, 0), (self.analog_fm_demod_cf_0, 0))
        self.connect((self.uhd_usrp_source_0, 0), (self.rational_resampler_xxx_0, 0))


    def closeEvent(self, event):
        self.settings = Qt.QSettings("GNU Radio", "manishkitry")
        self.settings.setValue("geometry", self.saveGeometry())
        self.stop()
        self.wait()

        event.accept()

    def get_volume(self):
        return self.volume

    def set_volume(self, volume):
        self.volume = volume

    def get_tuning(self):
        return self.tuning

    def set_tuning(self, tuning):
        self.tuning = tuning
        self.uhd_usrp_source_0.set_center_freq(self.tuning, 0)

    def get_sample_rate(self):
        return self.sample_rate

    def set_sample_rate(self, sample_rate):
        self.sample_rate = sample_rate
        self.uhd_usrp_source_0.set_samp_rate(self.sample_rate)

    def get_rf_decim(self):
        return self.rf_decim

    def set_rf_decim(self, rf_decim):
        self.rf_decim = rf_decim

    def get_interp(self):
        return self.interp

    def set_interp(self, interp):
        self.interp = interp

    def get_deviation(self):
        return self.deviation

    def set_deviation(self, deviation):
        self.deviation = deviation

    def get_audio_decim(self):
        return self.audio_decim

    def set_audio_decim(self, audio_decim):
        self.audio_decim = audio_decim

    def update_frequency(self):
        # time.sleep(5)
        current_freq = self.uhd_usrp_source_0.get_center_freq(0)
        new_freq = current_freq + 5e6  # Increment by 10 MHz
        self.uhd_usrp_source_0.set_center_freq(new_freq, 0)
        self.waterfall_display.set_center_freq(new_freq)
        print("Center frequency updated to:", new_freq)

    # def update_frequency(self):
    #     if self.toggle_frequency:
    #         new_freq = 98.3e6
    #     else:
    #         new_freq = 100.6e6

    #     self.uhd_usrp_source_0.set_center_freq(new_freq, 0)
    #     self.toggle_frequency = not self.toggle_frequency  # Toggle the flag

    #     print("Center frequency updated to:", new_freq)



def main(top_block_cls=manishkitry, options=None):

    qapp = Qt.QApplication(sys.argv)

    tb = top_block_cls()

    tb.start()

    tb.show()

    def sig_handler(sig=None, frame=None):
        tb.stop()
        tb.wait()

        Qt.QApplication.quit()

    signal.signal(signal.SIGINT, sig_handler)
    signal.signal(signal.SIGTERM, sig_handler)

    timer = Qt.QTimer()
    timer.start(500)
    timer.timeout.connect(lambda: None)

    qapp.exec_()

if __name__ == '__main__':
    main()
