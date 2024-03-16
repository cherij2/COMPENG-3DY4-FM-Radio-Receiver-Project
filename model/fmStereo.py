#
# Comp Eng 3DY4 (Computer Systems Integration Project)
#
# Copyright by Nicola Nicolici
# Department of Electrical and Computer Engineering
# McMaster University
# Ontario, Canada
#

import matplotlib.pyplot as plt
from scipy.io import wavfile
from scipy import signal
import numpy as np
import math

# use fmDemodArctan and fmPlotPSD
from fmSupportLib import fmDemodArctan, fmPlotPSD, derivativeDemod
from fmMonoBasic import LPF, BPF
# for take-home add your functions

rf_Fc = 100e3
rf_taps = 101


if_Fc = 16e3
if_taps = 101

audio_Fc = 16e3
audio_taps = 101
# add other settings for audio, like filter taps, ...

# flag that keeps track if your code is running for
# in-lab (il_vs_th = 0) vs takehome (il_vs_th = 1)
mode = 0
if (mode == 0):
	rf_Fs = 2.4e6
	rf_decim = 10
	if_Fs = 240e3
	audio_decim = 5
	audio_Fs = 48e3
elif(mode == 1):
	rf_Fs = 0.96e6
	rf_decim = 8
	if_Fs = 120e3
	audio_decim = 3
	audio_Fs = 40e3

il_vs_th = 0



def state_saving_convolution (h,xb, previous): #this function does convolution with state saving, replaces lfilter functionality

	len_xb = len(xb)
	len_h = len(h)
	yb = np.zeros(len_xb)
	for n in range(len(yb)):
		for k in range(len_h):
			if (n-k >= 0):
				yb[n] += h[k] * xb[n-k]
			else: #use the state from previous block
				yb[n] += h[k] * previous[n-k]
	previous  = xb[-len(h):-1]
	#now save the state for the next block
	
	return yb, previous



if __name__ == "__main__":

	# read the raw IQ data from the recorded file
	# IQ data is assumed to be in 8-bits unsigned (and interleaved)
	in_fname = "../data/samples0.raw"
	raw_data = np.fromfile(in_fname, dtype='uint8')
	print("Read raw RF data from \"" + in_fname + "\" in unsigned 8-bit format")
	# IQ data is normalized between -1 and +1 in 32-bit float format
	iq_data = (np.float32(raw_data) - 128.0)/128.0
	print("Reformatted raw RF data to 32-bit float format (" + str(iq_data.size * iq_data.itemsize) + " bytes)")

	# coefficients for the front-end low-pass filter
	rf_coeff = signal.firwin(rf_taps, rf_Fc/(rf_Fs/2), window=('hann'))

	# coefficients for the filter to extract mono audio
	if il_vs_th == 0:
		# to be updated by you during the in-lab session based on firwin
		# same principle  as for rf_coeff (but different arguments, of course)
		#copy pasted from lab1 blockProcessing
		audio_coeff = signal.firwin(if_taps, if_Fc/(if_Fs/2), window=('hann'))
		lpf_audio_coeff = LPF(if_Fs, if_Fc, if_taps)
		for x in range (len(audio_coeff)):
			#print("audio coeff with firwin: ", audio_coeff[x])
			print("own implementation audio coeff: ", lpf_audio_coeff[x])
		bpf_firwin_coeff = signal.firwin(if_taps, [18500/(if_Fs/2), 19500/(if_Fs/2)], pass_zero = False, fs = if_Fs )
		bpf_coeff = BPF(18500.0, 19500.0, if_Fs, if_taps)
		for x in range(len(bpf_coeff)):
			#print("BPF firwin coeffs: ", bpf_firwin_coeff[x])
			print("own implementation BPF: ", bpf_coeff[x])
		stereo_firwin_coeff = signal.firwin(if_taps, [22000/(if_Fs/2), 54000/(if_Fs/2)], pass_zero = False, fs = if_Fs )
		stereo_bpf_coeff = BPF(22000, 54000, if_Fs, if_taps)
		for x in range (len(stereo_bpf_coeff)):
			#print("Stereo band firwin coeffs: ", stereo_firwin_coeff[x])
			print("Stereo band own implementation: ", stereo_bpf_coeff[x])
	else:
		# to be updated by you for the takehome exercise
		# with your own code for impulse response generation
		audio_coeff = LPF(if_Fs, if_Fc, audio_taps)

	# set up the subfigures for plotting
	subfig_height = np.array([0.8, 2, 1.6]) # relative heights of the subfigures
	plt.rc('figure', figsize=(7.5, 7.5))	# the size of the entire figure
	fig, (ax0, ax1, ax2) = plt.subplots(nrows=3, gridspec_kw={'height_ratios': subfig_height})
	fig.subplots_adjust(hspace = .6)

	# select a block_size that is a multiple of KB
	# and a multiple of decimation factors
	block_size = 1024 * rf_decim * audio_decim * 2
	block_count = 0

	# states needed for continuity in block processing
	state_i_lpf_100k = np.zeros(rf_taps-1)
	state_q_lpf_100k = np.zeros(rf_taps-1)
	state_phase = 0
	I_prev = 0
	Q_prev = 0
	
	
	
	# add state as needed for the mono channel filter
	state_lpf = np.zeros(audio_taps-1)

	# audio buffer that stores all the audio blocks
	audio_data = np.array([]) # used to concatenate filtered blocks (audio data)

	# if the number of samples in the last block is less than the block size
	# it is fine to ignore the last few samples from the raw IQ file
	while (block_count+1)*block_size < 14:#len(iq_data):

		# if you wish to have shorter runtimes while troubleshooting
		# you can control the above loop exit condition as you see fit
		print('Processing block ' + str(block_count))

		#filter to extract the FM channel (I samples are even, Q samples are odd)
		i_filt, state_i_lpf_100k = signal.lfilter(rf_coeff, 1.0, \
				iq_data[(block_count)*block_size:(block_count+1)*block_size:2],
				zi=state_i_lpf_100k)
		q_filt, state_q_lpf_100k = signal.lfilter(rf_coeff, 1.0, \
				iq_data[(block_count)*block_size+1:(block_count+1)*block_size:2],
				zi=state_q_lpf_100k)


		# downsample the I/Q data from the FM channel
		i_ds = i_filt[::rf_decim]
		q_ds = q_filt[::rf_decim]

		# FM demodulator
		if il_vs_th == 0:
			# already given to you for the in-lab
			# take particular notice of the "special" state-saving
			fm_demod, state_phase = fmDemodArctan(i_ds, q_ds, state_phase)
			
		

		else:
			# you will need to implement your own FM demodulation based on:
			# https://www.embedded.com/dsp-tricks-frequency-demodulation-algorithms/
			# see more comments on fmSupportLib.py - take particular notice that
			# you MUST have also "custom" state-saving for your own FM demodulator
			fm_demod, state_phase,I_prev, Q_prev  = derivativeDemod(i_ds, q_ds, prev_phase = state_phase, I_prev=0.0, Q_prev=0.0 )






		# extract the mono audio data through filtering
		if il_vs_th == 0:
			# to be updated by you during the in-lab session based on lfilter
			# same principle as for i_filt or q_filt (but different arguments)
			audio_filt, state_lpf = signal.lfilter(audio_coeff, 1.0, fm_demod, zi = state_lpf)
		else:
		# 	# to be updated by you for the takehome exercise
		# 	# with your own code for BLOCK convolution
			audio_filt, state_lpf  = state_saving_convolution(audio_coeff,fm_demod, state_lpf)

		# downsample audio data
		# to be updated by you during in-lab (same code for takehome)
		audio_block = audio_filt[::audio_decim]

		# concatenate the most recently processed audio_block
		# to the previous blocks stored already in audio_data
		#
		audio_data = np.concatenate((audio_data, audio_block))
		#

		# to save runtime select the range of blocks to log data
		# this includes both saving binary files as well plotting PSD
		# below we assume we want to plot for graphs for blocks 10 and 11
		if block_count >= 10 and block_count < 12:

			# plot PSD of selected block after FM demodulation
			ax0.clear()
			fmPlotPSD(ax0, fm_demod, (rf_Fs/rf_decim)/1e3, subfig_height[0], \
					'Demodulated FM (block ' + str(block_count) + ')')
			# output binary file name (where samples are written from Python)
			fm_demod_fname = "../data/fm_demod_" + str(block_count) + ".bin"
			# create binary file where each sample is a 32-bit float
			fm_demod.astype('float32').tofile(fm_demod_fname)

			# plot PSD of selected block after extracting mono audio
			# ... change as needed
			fmPlotPSD(ax1, audio_filt, (rf_Fs/rf_decim)/1e3, subfig_height[1], 'Extracted Mono')

			# plot PSD of selected block after downsampling mono audio
			# ... change as needed
			fmPlotPSD(ax2, audio_block, audio_Fs/1e3, subfig_height[2], 'Downsampled Mono Audio')
			


			# save figure to file
			fig.savefig("../data/fmMonoBlock" + str(block_count) + ".png")

		block_count += 1
	
	print('Finished processing all the blocks from the recorded I/Q samples')

	# write audio data to file
	out_fname = "../data/fmMonoBlock.wav"
	wavfile.write(out_fname, int(audio_Fs), np.int16((audio_data/2)*32767))
	print("Written audio samples to \"" + out_fname + "\" in signed 16-bit format")

	# uncomment assuming you wish to show some plots
	#plt.show()