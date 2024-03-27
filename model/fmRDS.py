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
from fmPll import fmPll
from fmRRC import impulseResponseRootRaisedCosine
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

initial_state = {
	'integrator': 0.0,
    'phaseEst': 0.0,
    'feedbackI': 1.0,
    'feedbackQ': 0.0,
    'trigOffset': 0,
	'prev_ncoOut': 1.0
}

initial_state_rds = {
	'integrator': 0.0,
    'phaseEst': 0.0,
    'feedbackI': 1.0,
    'feedbackQ': 0.0,
    'trigOffset': 0,
	'prev_ncoOut': 1.0
}

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

def convolution_rs(yb, xb, h, ds, us, state):
    yb.clear()
    yb.extend([0.0] * ((len(xb) * us) // ds))

    for n in range(len(yb)):
        phase = (n * ds) % us
        for k in range(phase, len(h), us):
            dx = (ds * n - k) // us

            if dx >= 0:
                yb[n] += h[k] * xb[dx]
            else:
                yb[n] += h[k] * state[dx]

    state = xb[-len(state):]
    return yb, state

def delayBlock(input_block, state_block):
	output_block = np.concatenate((state_block, input_block[:-len(state_block)]))
	state_block = input_block[-len(state_block):]
	return output_block, state_block

if __name__ == "__main__":

	# read the raw IQ data from the recorded file
	# IQ data is assumed to be in 8-bits unsigned (and interleaved)
	in_fname = "../data/stereo_l0_r9.raw"
	raw_data = np.fromfile(in_fname, dtype='uint8')
	print("Read raw RF data from \"" + in_fname + "\" in unsigned 8-bit format")
	# IQ data is normalized between -1 and +1 in 32-bit float format
	iq_data = (np.float32(raw_data) - 128.0)/128.0
	print("Reformatted raw RF data to 32-bit float format (" + str(iq_data.size * iq_data.itemsize) + " bytes)")

	# coefficients for the front-end low-pass filter
	rf_coeff = signal.firwin(rf_taps, rf_Fc/(rf_Fs/2), window=('hann'))

	# coefficients for the filter to extract mono audio
	
	mono_coeff = signal.firwin(if_taps, if_Fc/(if_Fs/2), window=('hann'))

	pilot_coeffs = BPF(18500, 19500, if_Fs, if_taps)

	stereo_band_coeffs = BPF(22000, 54000, if_Fs, if_taps)
	
	rds_band_coeffs = BPF(54000, 60000, if_Fs, if_taps)

	rds_CR_coeffs = BPF(113500, 114500, if_Fs, if_taps)
    
	rds_DEM_LPF_coeffs = LPF(if_Fs, 3000, if_taps)
    
	rr_fs_out = 2375 * 30
	rr_fs_in = 240000
	rr_gcd = math.gcd(rr_fs_out, rr_fs_in)
	rr_us = rr_fs_out / rr_gcd
	rr_ds = rr_fs_in / rr_gcd
	rr_outp = [if_taps * rr_us]
	rr_fs = rr_us*if_Fs
	rr_fc = min((rr_us / rr_ds)*(rr_fs/2), rr_fs/2)
	rds_DEM_RR_coeffs = LPF(if_Fs, 3000, if_taps)

	crr_fs = (rr_us / rr_ds) * if_Fs

	# set up the subfigures for plotting
	# subfig_height = np.array([0.8, 2, 1.6]) # relative heights of the subfigures
	# plt.rc('figure', figsize=(7.5, 7.5))	# the size of the entire figure
	# fig, (ax0, ax1, ax2) = plt.subplots(nrows=3, gridspec_kw={'height_ratios': subfig_height})
	# fig.subplots_adjust(hspace = .6)

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
	
	
	state_delay = np.zeros((audio_taps-1)//2)
	rds_state_delay = np.zeros((audio_taps-1)//2)
	
	# add state as needed for the mono channel filter
	state_lpf = np.zeros(audio_taps-1)

	#state for stereoband and pilot 
	state_pilot = np.zeros(audio_taps-1)
	state_stereoband = np.zeros(audio_taps-1)
	bb_state_rds = np.zeros(audio_taps-1)
	CR_state_rds = np.zeros(audio_taps-1)
	CR_state_rds = np.zeros(audio_taps-1)
	CR_state_rds_pll = initial_state_rds
	rds_DEM_LPF_state = np.zeros(audio_taps-1)
	rds_DEM_rr_state = np.zeros(audio_taps-1)

	#mixer state
	state_mixer = np.zeros(audio_taps-1)
	# audio buffer that stores all the audio blocks
	mono_data = np.array([]) # used to concatenate filtered blocks (audio data)
	stereo_data = np.array([])
	left_data = np.array([])
	right_data = np.array([])
	state = initial_state

	# if the number of samples in the last block is less than the block size
	# it is fine to ignore the last few samples from the raw IQ file
	while (block_count+1)*block_size < 22*block_size: #len(iq_data):

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
		
		fm_demod, state_phase = fmDemodArctan(i_ds, q_ds, state_phase)
			
		# extract the mono audio data through filtering
		
		mono_filt, state_lpf = signal.lfilter(mono_coeff, 1.0, fm_demod, zi = state_lpf)
		
		final_mono, state_delay = delayBlock(mono_filt, state_delay)
		# downsample audio data
		# to be updated by you during in-lab (same code for takehome)fmPll		
		mono_block = final_mono[::audio_decim]

		# concatenate the most recently processed audio_block
		# to the previous blocks stored already in audio_data
		#
		mono_data = np.concatenate((mono_data, mono_block))
		#
		#convolution of fm demod and pilot coefficients
		pilot_filt, state_pilot = signal.lfilter(pilot_coeffs, 1.0, fm_demod, zi = state_pilot)

		#convolution of fm demod and stereband coefficients
		stereoband_filt, state_stereoband = signal.lfilter(stereo_band_coeffs, 1.0, fm_demod, zi = state_stereoband)


		ncoOut, state = fmPll(pilot_filt, 19000, if_Fs, 2.0, state = state)
		mixer =	2*ncoOut[:-1]*stereoband_filt

		mixer_filt, state_mixer = signal.lfilter(mono_coeff,1.0,mixer , zi = state_mixer)

		stereo_block = mixer_filt[::audio_decim]

		left_block = stereo_block + mono_block
		right_block = mono_block - stereo_block
		
		stereo_data = np.concatenate((stereo_data, stereo_block))


		left_data = np.concatenate((left_data, left_block))

		right_data = np.concatenate((right_data, right_block))

		#RDS IMPLEMENTATION

		#convolution for baseband
		bb_rds_filt, bb_state_rds = signal.lfilter(rds_band_coeffs, 1.0, fm_demod, zi = bb_state_rds) #CONVOLUTION BAND PASS

		delayed_rds, rds_state_delay = delayBlock(bb_rds_filt, rds_state_delay) #ALL PASS FILTER

		nonLin_rds = bb_rds_filt * bb_rds_filt # SQAURING NONLIN

		CR_rds_filt, CR_state_rds = signal.lfilter(rds_CR_coeffs, 1.0, nonLin_rds, zi = CR_state_rds) # BPF CARRIER RECOVERY

		rds_ncoOut, CR_state_rds_pll = fmPll(CR_rds_filt, 114000, if_Fs, ncoScale = 0.5, state = CR_state_rds_pll) # PLL

		#GRAPHING CARRIER RECVOERY 
		#t = np.arange(0, 1, 1/if_Fs)
		#t = np.arange(0, len(rds_ncoOut)/if_Fs, 1/if_Fs)
		if (block_count == 20):
			print(len(rds_ncoOut))
			plt.figure()
			plt.plot(rds_ncoOut[1850:2050], label='PLL Output')
			plt.plot(CR_rds_filt[1850:2050],label='PLL Input' )
			plt.xlabel('Time (s)')
			plt.ylabel('Amplitude')
			plt.legend(loc = 'upper left')
			plt.grid(True)
			plt.savefig("../data/fmPLLinput" + str(block_count) + ".png")
			plt.show()
			
		# mixer_rds = rds_ncoOut[:-1]*bb_rds_filt # MIX3R

    	# rds_DEM_LPF_filt, rds_DEM_LPF_state = signal.lfilter(rds_DEM_LPF_coeffs, 1.0, mixer_rds, zi = rds_DEM_LPF_state) #LPF
        
        # convolution_rs(rr_outp, rds_DEM_LPF_filt, rds_DEM_RR_coeffs, rr_ds, rr_us, rds_DEM_rr_state) #RATIONAL RESAMPLER

		# impulseResponseRootRaisedCosine(if_Fs, N_taps) # ROOT RAISED COSINE











  
		#print("mono filt size: ", len(mono_filt), " stereo filt size: ", len(mixer_filt))
		#print("mono block size: ", len(mono_block), " stereo block size: ", len(stereo_block))

		# to save runtime select the range of blocks to log data
		# this includes both saving binary files as well plotting PSD
		# below we assume we want to plot for graphs for blocks 10 and 11
		# if block_count >= 10 and block_count < 12:

		# 	# plot PSD of selected block after FM demodulation
		# 	ax0.clear()
		# 	fmPlotPSD(ax0, fm_demod, (rf_Fs/rf_decim)/1e3, subfig_height[0], \
		# 			'Demodulated FM (block ' + str(block_count) + ')')
		# 	# output binary file name (where samples are written from Python)
		# 	fm_demod_fname = "../data/fm_demod_" + str(block_count) + ".bin"
		# 	# create binary file where each sample is a 32-bit float
		# 	fm_demod.astype('float32').tofile(fm_demod_fname)

		# 	# plot PSD of selected block after extracting mono audio
		# 	# ... change as needed
		# 	fmPlotPSD(ax1, mixer_filt, (rf_Fs/rf_decim)/1e3, subfig_height[1], 'Extracted Stereo')

		# 	# plot PSD of selected block after downsampling mono audio
		# 	# ... change as needed
		# 	fmPlotPSD(ax2, stereo_block, audio_Fs/1e3, subfig_height[2], 'Downsampled Stereo Audio')
			


		# 	# # save figure to file
		# 	fig.savefig("../data/fmStereoBlock" + str(block_count) + ".png")

		block_count += 1
	
	# print('Finished processing all the blocks from the recorded I/Q samples')
	# print("length of right channel: ", len(right_data), " length of left channel", len(left_data))
	
	# audio_data = np.hstack((left_data[:, np.newaxis], right_data[:, np.newaxis]))
	# print("length of combined data: ", len(audio_data))
	# out_fname_left = "../data/fmleft.wav"
	# wavfile.write(out_fname_left, int(audio_Fs), np.int16((left_data/2)*32767))
	# print("Written audio samples to \"" + out_fname_left + "\" in signed 16-bit format")

	# out_fname_right = "../data/fmright.wav"
	# wavfile.write(out_fname_right, int(audio_Fs), np.int16((right_data/2)*32767))
	# print("Written audio samples to \"" + out_fname_right + "\" in signed 16-bit format")

	# # write audio data to file
	# out_fname = "../data/fmStereo.wav"
	# wavfile.write(out_fname, int(audio_Fs), np.int16((audio_data/2)*32767))
	# print("Written audio samples to \"" + out_fname + "\" in signed 16-bit format")

	# out_fname1 = "../data/fmMono.wav"
	# wavfile.write(out_fname1, int(audio_Fs), np.int16((mono_data/2)*32767))
	# print("Written audio samples to \"" + out_fname1 + "\" in signed 16-bit format")

	# uncomment assuming you wish to show some plots
	#plt.show()