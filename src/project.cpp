/*
Comp Eng 3DY4 (Computer Systems Integration Project)

Copyright by Nicola Nicolici
Department of Electrical and Computer Engineering
McMaster University
Ontario, Canada
*/

#include "dy4.h"
#include "filter.h"
#include "fourier.h"
#include "genfunc.h"
#include "iofunc.h"
#include "logfunc.h"
#include "RFfront.h"

// -------------------------- RF FRONT END ------------------------------

// function for computing the impulse response (reuse from previous experiment)


//creates blocks from input to pass into blockFIR, as stated to in the comments


// function to read audio data from a binary file that contains raw samples
// represented as 32-bit floats; we also assume two audio channels
// note: check the Python script that can prepare this type of files
// directly from .wav files


void readStdinBlockData(unsigned int num_samples, unsigned int block_id, std::vector<float> &block_data){
	std::vector<char> raw_data(num_samples);
	std::cin.read(reinterpret_cast<char*>(&raw_data[0]),num_samples*sizeof(char));
	for (int k=0; k<(int)num_samples; k++){
		block_data[k] = float((unsigned char)raw_data[k]-128/128.0);
	}
}
int main()
{
	
	/*
	# coefficients for the front-end low-pass filter
	rf_coeff = signal.firwin(rf_taps, rf_Fc/(rf_Fs/2), window=('hann'))

	# coefficients for the filter to extract mono audio
	if il_vs_th == 0:
		# to be updated by you during the in-lab session based on firwin
		# same principle  as for rf_coeff (but different arguments, of course)
		#copy pasted from lab1 blockProcessing
		audio_coeff = signal.firwin(audio_taps, audio_Fc/(audio_Fs/2), window=('hann'))
	else:
		# to be updated by you for the takehome exercise
		# with your own code for impulse response generation
		audio_coeff = LPF(audio_Fs, audio_Fc, audio_taps)

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
	
	 */

	
	float RF_Fs = 2400e3;
	float RF_Fc = 100e3;
	float IF_Fs = 240e3;
	float mono_Fc = 16e3;
	float num_Taps = 101;
	int rf_decim = 10;
	int audio_decim = 5;
	std::vector<float> RF_h;
	std::vector<float> IF_h;
	
	std::vector<float> i_data, q_data;
	std::vector<float> filt_i, filt_q;
	std::vector<float> demod;
	float prev_i = 0.0; 
	float prev_q =0.0;
	std::vector<float> mono_out;


	
	int BLOCK_SIZE = 1024 * rf_decim * audio_decim * 2;

	for (unsigned int block_id = 0;  ; block_id++) {
		std::vector<float> block_data(BLOCK_SIZE);
		readStdinBlockData(BLOCK_SIZE, block_id, block_data);
		if((std::cin.rdstate()) != 0) {
			std::cerr << "End of input stream reached" << std::endl;
			exit(1);
		}
		std::cerr << "Read block " << block_id << std::endl;
		split_audio_iq(block_data, i_data, q_data);
		impulseResponseLPF(RF_Fs, RF_Fc, num_Taps, RF_h);
		conv_ds_slow(filt_i, i_data, RF_h, rf_decim);
		conv_ds_slow(filt_q, q_data, RF_h, rf_decim);
		FM_demod(filt_i, filt_q, prev_i, prev_q, demod);
		impulseResponseLPF(IF_Fs, mono_Fc, num_Taps, IF_h);
		conv_ds_slow(mono_out, demod, IF_h, audio_decim);

	}

	

	return 0;
}
