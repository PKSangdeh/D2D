#!/usr/bin/env python2
# -*- coding: utf-8 -*-
##################################################
# GNU Radio Python Flow Graph
# Title: Top Block
# Generated: Tue Jul 16 22:14:02 2019
##################################################

if __name__ == '__main__':
    import ctypes
    import sys
    if sys.platform.startswith('linux'):
        try:
            x11 = ctypes.cdll.LoadLibrary('libX11.so')
            x11.XInitThreads()
        except:
            print "Warning: failed to XInitThreads()"

from PyQt4 import Qt
from gnuradio import analog
from gnuradio import blocks
from gnuradio import eng_notation
from gnuradio import gr
from gnuradio import qtgui
from gnuradio import uhd
from gnuradio.eng_option import eng_option
from gnuradio.filter import firdes
from gnuradio.qtgui import Range, RangeWidget
from optparse import OptionParser
import sip
import sys
import time


class top_block(gr.top_block, Qt.QWidget):

    def __init__(self):
        gr.top_block.__init__(self, "Top Block")
        Qt.QWidget.__init__(self)
        self.setWindowTitle("Top Block")
        try:
            self.setWindowIcon(Qt.QIcon.fromTheme('gnuradio-grc'))
        except:
            pass
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

        self.settings = Qt.QSettings("GNU Radio", "top_block")
        self.restoreGeometry(self.settings.value("geometry").toByteArray())

        ##################################################
        # Variables
        ##################################################
        self.samp_rate = samp_rate = 5000000
        self.rf_gain = rf_gain = 10
        self.phase_offset3 = phase_offset3 = 0
        self.phase_offset2 = phase_offset2 = 0
        self.phase_offset1 = phase_offset1 = 0
        self.fft_size = fft_size = 4096
        self.const = const = 0.2
        self.cent_freq = cent_freq = 2100000000
        self.bandwidth = bandwidth = samp_rate

        ##################################################
        # Blocks
        ##################################################
        self._phase_offset3_range = Range(-1.5, 1.5, 0.02, 0, 200)
        self._phase_offset3_win = RangeWidget(self._phase_offset3_range, self.set_phase_offset3, "phase_offset3", "counter_slider", float)
        self.top_layout.addWidget(self._phase_offset3_win)
        self._phase_offset2_range = Range(-1.5, 1.5, 0.02, 0, 200)
        self._phase_offset2_win = RangeWidget(self._phase_offset2_range, self.set_phase_offset2, "phase_offset2", "counter_slider", float)
        self.top_layout.addWidget(self._phase_offset2_win)
        self._phase_offset1_range = Range(-1.5, 1.5, 0.02, 0, 200)
        self._phase_offset1_win = RangeWidget(self._phase_offset1_range, self.set_phase_offset1, "phase_offset1", "counter_slider", float)
        self.top_layout.addWidget(self._phase_offset1_win)
        self.uhd_usrp_source_0_0 = uhd.usrp_source(
        	",".join(("addr=192.168.1.185", "")),
        	uhd.stream_args(
        		cpu_format="fc32",
        		channels=range(1),
        	),
        )
        self.uhd_usrp_source_0_0.set_samp_rate(samp_rate)
        self.uhd_usrp_source_0_0.set_center_freq(cent_freq, 0)
        self.uhd_usrp_source_0_0.set_gain(28, 0)
        self.uhd_usrp_source_0_0.set_antenna("TX/RX", 0)
        self.uhd_usrp_source_0_0.set_bandwidth(bandwidth, 0)
        self.uhd_usrp_sink_0_0 = uhd.usrp_sink(
        	",".join(("addr0=192.168.1.175, addr1=192.168.1.176,addr2=192.168.1.177,addr3=192.168.1.178", "")),
        	uhd.stream_args(
        		cpu_format="fc32",
        		channels=range(4),
        	),
        )
        self.uhd_usrp_sink_0_0.set_clock_source("external", 0)
        self.uhd_usrp_sink_0_0.set_time_source("external", 0)
        self.uhd_usrp_sink_0_0.set_clock_source("external", 1)
        self.uhd_usrp_sink_0_0.set_time_source("external", 1)
        self.uhd_usrp_sink_0_0.set_clock_source("external", 2)
        self.uhd_usrp_sink_0_0.set_time_source("external", 2)
        self.uhd_usrp_sink_0_0.set_clock_source("external", 3)
        self.uhd_usrp_sink_0_0.set_time_source("external", 3)
        self.uhd_usrp_sink_0_0.set_time_unknown_pps(uhd.time_spec())
        self.uhd_usrp_sink_0_0.set_samp_rate(samp_rate)
        self.uhd_usrp_sink_0_0.set_center_freq(cent_freq, 0)
        self.uhd_usrp_sink_0_0.set_gain(rf_gain, 0)
        self.uhd_usrp_sink_0_0.set_antenna("TX/RX", 0)
        self.uhd_usrp_sink_0_0.set_bandwidth(bandwidth, 0)
        self.uhd_usrp_sink_0_0.set_center_freq(cent_freq, 1)
        self.uhd_usrp_sink_0_0.set_gain(rf_gain, 1)
        self.uhd_usrp_sink_0_0.set_antenna("TX/RX", 1)
        self.uhd_usrp_sink_0_0.set_bandwidth(bandwidth, 1)
        self.uhd_usrp_sink_0_0.set_center_freq(cent_freq, 2)
        self.uhd_usrp_sink_0_0.set_gain(rf_gain, 2)
        self.uhd_usrp_sink_0_0.set_antenna("TX/RX", 2)
        self.uhd_usrp_sink_0_0.set_bandwidth(bandwidth, 2)
        self.uhd_usrp_sink_0_0.set_center_freq(cent_freq, 3)
        self.uhd_usrp_sink_0_0.set_gain(rf_gain, 3)
        self.uhd_usrp_sink_0_0.set_bandwidth(bandwidth, 3)
        self.qtgui_sink_x_0_0_0 = qtgui.sink_c(
        	fft_size, #fftsize
        	firdes.WIN_BLACKMAN_hARRIS, #wintype
        	cent_freq, #fc
        	samp_rate, #bw
        	"177", #name
        	True, #plotfreq
        	True, #plotwaterfall
        	True, #plottime
        	True, #plotconst
        )
        self.qtgui_sink_x_0_0_0.set_update_time(1.0/10)
        self._qtgui_sink_x_0_0_0_win = sip.wrapinstance(self.qtgui_sink_x_0_0_0.pyqwidget(), Qt.QWidget)
        self.top_layout.addWidget(self._qtgui_sink_x_0_0_0_win)
        
        self.qtgui_sink_x_0_0_0.enable_rf_freq(False)
        
        
          
        self.blocks_skiphead_0_0_0 = blocks.skiphead(gr.sizeof_gr_complex*1, 1*samp_rate)
        self.blocks_multiply_xx_0_0_0_1 = blocks.multiply_vcc(1)
        self.blocks_multiply_xx_0_0_0_0 = blocks.multiply_vcc(1)
        self.blocks_multiply_xx_0_0_0 = blocks.multiply_vcc(1)
        self.blocks_multiply_const_vxx_0_0_0_0_0_1 = blocks.multiply_const_vcc((const, ))
        self.blocks_multiply_const_vxx_0_0_0_0_0_0 = blocks.multiply_const_vcc((const, ))
        self.blocks_multiply_const_vxx_0_0_0_0_0 = blocks.multiply_const_vcc((const, ))
        self.blocks_multiply_const_vxx_0_0 = blocks.multiply_const_vcc((const, ))
        self.blocks_magphase_to_complex_0_0_0_1 = blocks.magphase_to_complex(1)
        self.blocks_magphase_to_complex_0_0_0_0 = blocks.magphase_to_complex(1)
        self.blocks_magphase_to_complex_0_0_0 = blocks.magphase_to_complex(1)
        self.blocks_file_source_0_1_0_0_0_0_0 = blocks.file_source(gr.sizeof_gr_complex*1, "/home/sdr/Desktop/Pedram/D2D_wo_ext_clk/signals/nr_bf_tx_3.dat", True)
        self.blocks_file_source_0_1_0_0_0_0 = blocks.file_source(gr.sizeof_gr_complex*1, "/home/sdr/Desktop/Pedram/D2D_wo_ext_clk/signals/nr_bf_tx_2.dat", True)
        self.blocks_file_source_0_1_0_0_0 = blocks.file_source(gr.sizeof_gr_complex*1, "/home/sdr/Desktop/Pedram/D2D_wo_ext_clk/signals/nr_bf_tx_1.dat", True)
        self.blocks_file_source_0_1_0_0 = blocks.file_source(gr.sizeof_gr_complex*1, "/home/sdr/Desktop/Pedram/D2D_wo_ext_clk/signals/nr_bf_tx_0.dat", True)
        self.analog_const_source_x_0_1_0_1 = analog.sig_source_f(0, analog.GR_CONST_WAVE, 0, 0, 1)
        self.analog_const_source_x_0_1_0_0 = analog.sig_source_f(0, analog.GR_CONST_WAVE, 0, 0, 1)
        self.analog_const_source_x_0_1_0 = analog.sig_source_f(0, analog.GR_CONST_WAVE, 0, 0, 1)
        self.analog_const_source_x_0_0_0_0_1 = analog.sig_source_f(0, analog.GR_CONST_WAVE, 0, 0, phase_offset1)
        self.analog_const_source_x_0_0_0_0_0 = analog.sig_source_f(0, analog.GR_CONST_WAVE, 0, 0, phase_offset2)
        self.analog_const_source_x_0_0_0_0 = analog.sig_source_f(0, analog.GR_CONST_WAVE, 0, 0, phase_offset3)

        ##################################################
        # Connections
        ##################################################
        self.connect((self.analog_const_source_x_0_0_0_0, 0), (self.blocks_magphase_to_complex_0_0_0, 1))    
        self.connect((self.analog_const_source_x_0_0_0_0_0, 0), (self.blocks_magphase_to_complex_0_0_0_0, 1))    
        self.connect((self.analog_const_source_x_0_0_0_0_1, 0), (self.blocks_magphase_to_complex_0_0_0_1, 1))    
        self.connect((self.analog_const_source_x_0_1_0, 0), (self.blocks_magphase_to_complex_0_0_0, 0))    
        self.connect((self.analog_const_source_x_0_1_0_0, 0), (self.blocks_magphase_to_complex_0_0_0_0, 0))    
        self.connect((self.analog_const_source_x_0_1_0_1, 0), (self.blocks_magphase_to_complex_0_0_0_1, 0))    
        self.connect((self.blocks_file_source_0_1_0_0, 0), (self.blocks_multiply_const_vxx_0_0, 0))    
        self.connect((self.blocks_file_source_0_1_0_0_0, 0), (self.blocks_multiply_xx_0_0_0_1, 0))    
        self.connect((self.blocks_file_source_0_1_0_0_0_0, 0), (self.blocks_multiply_xx_0_0_0_0, 0))    
        self.connect((self.blocks_file_source_0_1_0_0_0_0_0, 0), (self.blocks_multiply_xx_0_0_0, 0))    
        self.connect((self.blocks_magphase_to_complex_0_0_0, 0), (self.blocks_multiply_xx_0_0_0, 1))    
        self.connect((self.blocks_magphase_to_complex_0_0_0_0, 0), (self.blocks_multiply_xx_0_0_0_0, 1))    
        self.connect((self.blocks_magphase_to_complex_0_0_0_1, 0), (self.blocks_multiply_xx_0_0_0_1, 1))    
        self.connect((self.blocks_multiply_const_vxx_0_0, 0), (self.uhd_usrp_sink_0_0, 0))    
        self.connect((self.blocks_multiply_const_vxx_0_0_0_0_0, 0), (self.uhd_usrp_sink_0_0, 3))    
        self.connect((self.blocks_multiply_const_vxx_0_0_0_0_0_0, 0), (self.uhd_usrp_sink_0_0, 2))    
        self.connect((self.blocks_multiply_const_vxx_0_0_0_0_0_1, 0), (self.uhd_usrp_sink_0_0, 1))    
        self.connect((self.blocks_multiply_xx_0_0_0, 0), (self.blocks_multiply_const_vxx_0_0_0_0_0, 0))    
        self.connect((self.blocks_multiply_xx_0_0_0_0, 0), (self.blocks_multiply_const_vxx_0_0_0_0_0_0, 0))    
        self.connect((self.blocks_multiply_xx_0_0_0_1, 0), (self.blocks_multiply_const_vxx_0_0_0_0_0_1, 0))    
        self.connect((self.blocks_skiphead_0_0_0, 0), (self.qtgui_sink_x_0_0_0, 0))    
        self.connect((self.uhd_usrp_source_0_0, 0), (self.blocks_skiphead_0_0_0, 0))    

    def closeEvent(self, event):
        self.settings = Qt.QSettings("GNU Radio", "top_block")
        self.settings.setValue("geometry", self.saveGeometry())
        event.accept()


    def get_samp_rate(self):
        return self.samp_rate

    def set_samp_rate(self, samp_rate):
        self.samp_rate = samp_rate
        self.set_bandwidth(self.samp_rate)
        self.uhd_usrp_source_0_0.set_samp_rate(self.samp_rate)
        self.qtgui_sink_x_0_0_0.set_frequency_range(self.cent_freq, self.samp_rate)
        self.uhd_usrp_sink_0_0.set_samp_rate(self.samp_rate)

    def get_rf_gain(self):
        return self.rf_gain

    def set_rf_gain(self, rf_gain):
        self.rf_gain = rf_gain
        self.uhd_usrp_source_0_0.set_gain(self.rf_gain, 1)
        	
        self.uhd_usrp_source_0_0.set_gain(self.rf_gain, 2)
        	
        self.uhd_usrp_source_0_0.set_gain(self.rf_gain, 3)
        	
        self.uhd_usrp_sink_0_0.set_gain(self.rf_gain, 0)
        	
        self.uhd_usrp_sink_0_0.set_gain(self.rf_gain, 1)
        	
        self.uhd_usrp_sink_0_0.set_gain(self.rf_gain, 2)
        	
        self.uhd_usrp_sink_0_0.set_gain(self.rf_gain, 3)
        	

    def get_phase_offset3(self):
        return self.phase_offset3

    def set_phase_offset3(self, phase_offset3):
        self.phase_offset3 = phase_offset3
        self.analog_const_source_x_0_0_0_0.set_offset(self.phase_offset3)

    def get_phase_offset2(self):
        return self.phase_offset2

    def set_phase_offset2(self, phase_offset2):
        self.phase_offset2 = phase_offset2
        self.analog_const_source_x_0_0_0_0_0.set_offset(self.phase_offset2)

    def get_phase_offset1(self):
        return self.phase_offset1

    def set_phase_offset1(self, phase_offset1):
        self.phase_offset1 = phase_offset1
        self.analog_const_source_x_0_0_0_0_1.set_offset(self.phase_offset1)

    def get_fft_size(self):
        return self.fft_size

    def set_fft_size(self, fft_size):
        self.fft_size = fft_size

    def get_const(self):
        return self.const

    def set_const(self, const):
        self.const = const
        self.blocks_multiply_const_vxx_0_0.set_k((self.const, ))
        self.blocks_multiply_const_vxx_0_0_0_0_0.set_k((self.const, ))
        self.blocks_multiply_const_vxx_0_0_0_0_0_0.set_k((self.const, ))
        self.blocks_multiply_const_vxx_0_0_0_0_0_1.set_k((self.const, ))

    def get_cent_freq(self):
        return self.cent_freq

    def set_cent_freq(self, cent_freq):
        self.cent_freq = cent_freq
        self.uhd_usrp_source_0_0.set_center_freq(self.cent_freq, 0)
        self.uhd_usrp_source_0_0.set_center_freq(self.cent_freq, 1)
        self.uhd_usrp_source_0_0.set_center_freq(self.cent_freq, 2)
        self.uhd_usrp_source_0_0.set_center_freq(self.cent_freq, 3)
        self.qtgui_sink_x_0_0_0.set_frequency_range(self.cent_freq, self.samp_rate)
        self.uhd_usrp_sink_0_0.set_center_freq(self.cent_freq, 0)
        self.uhd_usrp_sink_0_0.set_center_freq(self.cent_freq, 1)
        self.uhd_usrp_sink_0_0.set_center_freq(self.cent_freq, 2)
        self.uhd_usrp_sink_0_0.set_center_freq(self.cent_freq, 3)

    def get_bandwidth(self):
        return self.bandwidth

    def set_bandwidth(self, bandwidth):
        self.bandwidth = bandwidth
        self.uhd_usrp_source_0_0.set_bandwidth(self.bandwidth, 0)
        self.uhd_usrp_source_0_0.set_bandwidth(self.bandwidth, 1)
        self.uhd_usrp_source_0_0.set_bandwidth(self.bandwidth, 2)
        self.uhd_usrp_source_0_0.set_bandwidth(self.bandwidth, 3)
        self.uhd_usrp_sink_0_0.set_bandwidth(self.bandwidth, 0)
        self.uhd_usrp_sink_0_0.set_bandwidth(self.bandwidth, 1)
        self.uhd_usrp_sink_0_0.set_bandwidth(self.bandwidth, 2)
        self.uhd_usrp_sink_0_0.set_bandwidth(self.bandwidth, 3)


def main(top_block_cls=top_block, options=None):

    from distutils.version import StrictVersion
    if StrictVersion(Qt.qVersion()) >= StrictVersion("4.5.0"):
        style = gr.prefs().get_string('qtgui', 'style', 'raster')
        Qt.QApplication.setGraphicsSystem(style)
    qapp = Qt.QApplication(sys.argv)

    tb = top_block_cls()
    tb.start()
    tb.show()

    def quitting():
        tb.stop()
        tb.wait()
    qapp.connect(qapp, Qt.SIGNAL("aboutToQuit()"), quitting)
    qapp.exec_()


if __name__ == '__main__':
    main()
