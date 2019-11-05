# D2D
Description:

   - This codes enables the D2D communication at cellular networks. The network has two subsystems, cellular network and D2D network.
     The cellular network includes a multi-antenna BS and multiple single-antenna CUEs. The D2D network includes two DUEs, each of 
     which may have one or multiple-antennas. But, make sure there is enough DoF to cancell cross-network interference. For example,
     a BS with 4 antennas, 2 single-antenna CUEs, DUE1 with 2 antennas, and DUE2 with one antenna is a proper configuration. The D2D
     and cellular network work simultanusly in the same band and the ultimate objective of this prohect is enhancing the spectral 
     efficiency.


PHY Type:

   - The Cellular network always uses 5G NR-like PHY with numerology mu=0 with 512 fft points, 15kHz subcarrier spacing, normal CP, 
     and 1ms/slot duration.
   - The D2D network can employ either 5G NR-like PHY with aforementioned parameters or IEEE 802.11 PHY with 64 subcarriers, 8 CP
     length, 4 pilot subcarriers, and 48 payload subcarriers.


Procedure at high level:

   + Preliminary phase: BS and DUE1 performs calibration to calculate uplink/downlink channel mismatch and maintain reciprosity.

   + Uplink phase:   CUEs and DUE2 send their streams and DUE1 and BS recover the desired streams and design spatial filters to 
                     pre-cancel cross-network interference in the downlink phase.
   + Downlink phase: During downlink phase, DUE1 and BS pre-code the data with the spatial filters obtained in the previous phase.
                     Ideally, CUEs and DUE2 should not receive any interference and decode their data in regular manner.


Functionality of some of the codes:

   + Any code with _encode.m at the end will generate proper files to be transmitted for calibratio or data transmission
   + All codes having calib in their name are a part of calibration procedure.
   + PN_generator creates a preamble for time synchornization. It can be zadoff-chu or pseud-noise sequence.
   + Pattern_generator defines different DMRS pattern. It has many pre-defined pattern and also you can add your desired pattern.

Purpose of stored .mat files:

   + pn stores the preamble
   + G_stack has the spatial filters for each subcarrier.
   + Patterns stores all the user-defined patterns for DMRS.
   + D_calib and C_calib store the calibration coefficients for DUE1 and BS, respectively.


Order of offine steps:

   + Calibration ----> Uplink ----> Uplink decoding ----> BF design ----> Precoding ----> Downlink Tx ----> User decodes their signals

Outputs of different codes:

   + Based on function of each code, a code may return EVM, synchornization result, freq. offset, constellation of decoded signal, 
     calibration coefficient, etc.
