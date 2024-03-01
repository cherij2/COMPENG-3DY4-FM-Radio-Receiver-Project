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

// -------------------------- RF FRONT END ------------------------------

// function for computing the impulse response (reuse from previous experiment)
void LPF(float Fs, float Fc, unsigned short int num_taps, std::vector<float> &h)
{
	// allocate memory for the impulse response
	h.clear();
	h.resize(num_taps, 0.0);
	float Norm_cutoff = Fc / (Fs / 2);

	for (int m = 0; m < num_taps; m++)
	{
		if (m == ((num_taps - 1) / 2))
		{
			h[m] = Norm_cutoff;
		}
		else
		{
			h[m] = Norm_cutoff * (std::sin(PI * Norm_cutoff * (m - (num_taps - 1) / 2)) / (PI * Norm_cutoff * (m - (num_taps - 1) / 2)));
		}
		h[m] = h[m] * pow((std::sin((m * PI) / num_taps)), 2);
	}
}


//creates blocks from input to pass into blockFIR, as stated to in the comments
void createBlocks(const std::vector<float> &x, int block_size, std::vector<std::vector<float>> &blocks)
{
	blocks.clear();
	int num_blocks = x.size() / block_size;//calculate number of blocks based on the size of input vector x
	for (int i = 0; i < num_blocks; i++)//iterate over number of blocks
	{
		std::vector<float> block(x.begin() + i * block_size, x.begin() + (i + 1) * block_size);
		blocks.push_back(block);//add the current block to blocks vector
	}//divides input vector into block_size number of blocks, this is done for processing data in blocks
}


void blockConvolutionFIR(std::vector<float> & yb, const std::vector<float> &xb, const std::vector<float> &h, std::vector<float> &state)
{//parameters include yb which is output block, xb input signal, h is the impulse response of the filter, state is state that will be used and updated to be used for the next block
	yb.clear();//this implementation copies the python code
	yb.resize(xb.size(), 0.0);

	for (int n = 0; n < xb.size(); n++)
	{
		for (int k = 0; k < h.size(); k++)
		{
			if (n - k >= 0)
			{
				yb[n] += h[k] * xb[n - k];
			}
			else
			{
				yb[n] += h[k] * state[h.size() - 1 + (n - k)];//
			}
		}
	}

	for (int i = 0; i < h.size() - 1; i++)//updates the state at the end for the next block to be used
	{
		state[i] = xb[xb.size() - h.size() + 1 + i];
	}
}

// function to read audio data from a binary file that contains raw samples
// represented as 32-bit floats; we also assume two audio channels
// note: check the Python script that can prepare this type of files
// directly from .wav files
void read_raw_data(const std::string in_fname, std::vector<float> &raw_data)
{
	// file descriptor for the input to be read
	std::ifstream fdin(in_fname, std::ios::binary);
	if (!fdin)
	{
		std::cout << "File " << in_fname << " not found ... exiting\n";
		exit(1);
	}
	else
	{
		std::cout << "Reading raw audio from \"" << in_fname << "\"\n";
	}
	// search for end of file to count the number of samples to be read
	fdin.seekg(0, std::ios::end);
	// we assume the Python script has written data in 32-bit floats
	const unsigned int num_samples = fdin.tellg() / sizeof(float);

	// allocate memory space to store all the samples
	raw_data.clear();
	raw_data.resize(num_samples);
	// back to the beginning of the file to read all samples at once
	fdin.seekg(0, std::ios::beg);
	// do a single read for audio data from the input file stream
	fdin.read(reinterpret_cast<char *>(&raw_data[0]),
				num_samples * sizeof(float));
	// close the input file
	fdin.close();
}

// function to split an audio data where the left channel is in even samples
// and the right channel is in odd samples
// (I samples are even, Q samples are odd)
void split_raw_data_into_channels(const std::vector<float> &raw_data, std::vector<float> &i_data, std::vector<float> &q_data)
{
	for (int i = 0; i < (int)raw_data.size(); i++)
	{
		if (i % 2 == 0)
			i_data.push_back(raw_data[i]);
		else
			q_data.push_back(raw_data[i]);
	}
}

// function to write audio data to a binary file that contains raw samples
// represented as 32-bit floats; we also assume two audio channels
// note: check the python script that can read this type of files
// and then reformat them to .wav files to be run on third-party players
void write_raw_data(const std::string out_fname, const std::vector<float> &i_data, const std::vector<float> &q_data)
{
	// file descriptor for the output to be written
	if (i_data.size() != q_data.size())
	{
		std::cout << "Something got messed up with audio channels\n";
		std::cout << "They must have the same size ... exiting\n";
		exit(1);
	}
	else
	{
		std::cout << "Writing raw audio to \"" << out_fname << "\"\n";
	}
	std::ofstream fdout(out_fname, std::ios::binary);
	for (int i = 0; i < (int)i_data.size(); i++)
	{
		// we assume we have handled a stereo audio file
		// hence, we must interleave the two channels
		// (change as needed if testing with mono files)
		fdout.write(reinterpret_cast<const char *>(&i_data[i]),
					sizeof(i_data[i]));
		fdout.write(reinterpret_cast<const char *>(&q_data[i]),
					sizeof(q_data[i]));
	}
	fdout.close();
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
	
	float RF_M0_Fs = 2400;
	float RF_M0_Fc = 16,000;
	float RF_M0_numTaps = 101;
	std::vector<float> RF_M0_h;
	
	std::vector<float> raw_data;
	std::vector<float> i_data(), q_data();
	std::vector<float> iq_data(raw_data.size());
	
	const std::string in_fname = "../data/samples3.raw";
	std::vector<float> raw_data;
	read_raw_data(in_fname, raw_data);
	
	for (int i = 0; i < raw_data.size(); ++i) {
        iq_data[i] = (static_cast<float>(raw_data[i]) - 128.0) / 128.0;
    }
    
    rf_coeff = LPF(RF_M0_Fs, RF_M0_Fc, RF_M0_numTaps, RF_M0_h);
    
	
	
	
	
	// binary files can be generated through the
	// Python models from the "../model/" sub-folder
	const std::string in_fname = "../data/fm_demod_10.bin";
	std::vector<float> bin_data;
	readBinData(in_fname, bin_data);

	// generate an index vector to be used by logVector on the X axis
	std::vector<float> vector_index;
	genIndexVector(vector_index, bin_data.size());
	// log time data in the "../data/" subfolder in a file with the following name
	// note: .dat suffix will be added to the log file in the logVector function
	logVector("demod_time", vector_index, bin_data);

	// take a slice of data with a limited number of samples for the Fourier transform
	// note: NFFT constant is actually just the number of points for the
	// Fourier transform - there is no FFT implementation ... yet
	// unless you wish to wait for a very long time, keep NFFT at 1024 or below
	std::vector<float> slice_data = \
		std::vector<float>(bin_data.begin(), bin_data.begin() + NFFT);
	// note: make sure that binary data vector is big enough to take the slice

	// declare a vector of complex values for DFT
	std::vector<std::complex<float>> Xf;
	// ... in-lab ...
	// compute the Fourier transform
	// the function is already provided in fourier.cpp

	// compute the magnitude of each frequency bin
	// note: we are concerned only with the magnitude of the frequency bin
	// (there is NO logging of the phase response)
	std::vector<float> Xmag;
	// ... in-lab ...
	// compute the magnitude of each frequency bin
	// the function is already provided in fourier.cpp

	// log the frequency magnitude vector
	vector_index.clear();
	genIndexVector(vector_index, Xmag.size());
	logVector("demod_freq", vector_index, Xmag); // log only positive freq

	// for your take-home exercise - repeat the above after implementing
	// your OWN function for PSD based on the Python code that has been provided
	// note the estimate PSD function should use the entire block of "bin_data"
	//
	// ... complete as part of the take-home ...
	//

	// if you wish to write some binary files, see below example
	//
	// const std::string out_fname = "../data/outdata.bin";
	// writeBinData(out_fname, bin_data);
	//
	// output files can be imported, for example, in Python
	// for additional analysis or alternative forms of visualization
	
	
	
	
	
	// naturally, you can comment the line below once you are comfortable to run GNU plot
	std::cout << "Run: gnuplot -e 'set terminal png size 1024,768' ../data/example.gnuplot > ../data/example.png\n";

	return 0;
}
