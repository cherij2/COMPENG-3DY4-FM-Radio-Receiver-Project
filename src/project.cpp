

// -------------------------- RF FRONT END ------------------------------




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



void readStdinBlockData(unsigned int num_samples, unsigned int block_id, std::vector<float> &block_data){
	std::vector<char> raw_data(num_samples);
	std::cin.read(reinterpret_cast<char*>(&raw_data[0]), num_samples*sizeof(char));
	for (int k=0; k<(int)num_samples; k++){
		block_data[k] = float(((unsigned char)raw_data[k]-128)/128.0);
	}
}

/* int main(int argc, char* argv[])
{
	// Default mode 0
	int mode = 0;
	std::string audio_channel = "m"; // default to mono

	// Mode Selection
	if (argc < 2) {
		std::cerr << "Operating in default mode 0 with mono audio" << std::endl;
	} else if (argc == 2) {
		mode = atoi(argv[1]);
		if (mode > 3) {
			std::cerr << "Wrong mode " << mode << ". Mode must be between 0 and 3." << std::endl;
			exit(1);
		}
		std::cerr << "Operating in mode " << mode << " with mono audio" << std::endl;
	} else if (argc == 3) {
		mode = atoi(argv[1]);
		audio_channel = argv[2];

		if (mode > 3 || (audio_channel != "m" && audio_channel != "s" && audio_channel != "r")) {
			std::cerr << "Usage: " << argv[0] << " <mode> [audio channel]" << std::endl;
			std::cerr << "<mode> is a value from 0-3 and [audio channel] can be 'm' for mono, 's' for stereo, or 'r' for RDS." << std::endl;
			exit(1);
		}
		std::cerr << "Operating in mode " << mode << " with audio channel " << audio_channel << std::endl;
	} else {
		std::cerr << "Usage: " << argv[0] << " <mode> [audio channel]" << std::endl;
		std::cerr << "<mode> is a value from 0-3 and [audio channel] can be 'm' for mono, 's' for stereo, or 'r' for RDS." << std::endl;
		exit(1);
	} *///above implements the 3 arguments that will be given ./project <mode> [audio_channel]

int main(int argc, char *argv[])
{

	// Default mode 0
	int mode = 0;

	// Mode Selection
	if (argc<2){
		std::cerr << "Operating in default mode 0" << std::endl;
	} else if (argc==2){
		mode=atoi(argv[1]);
		if (mode>3){
			std::cerr << "Wrong mode " << mode << std::endl;
			exit(1);
		}
	} else {
		std::cerr << "Usage: " << argv[0] << std::endl;
		std::cerr << "or " << std::endl;
		std::cerr << "Usage: " << argv[0] << " <mode>" << std::endl;
		std::cerr << "\t\t <mode> is a value from 0-3" << argv[0] << std::endl;
		exit(1);
	}
	std::cerr << "Operating in mode " << mode << std::endl;


	float RF_Fs = 2400e3;
	float RF_Fc = 100e3;
	float IF_Fs = 240e3;
	float mono_Fc = 16e3;
	float num_Taps = 151;
	int rf_decim = 10;
	int audio_decim = 5;
	int audio_expan;
	int BLOCK_SIZE = 1024*rf_decim*audio_decim*2;
	if (mode == 1){
		RF_Fs = 960e3;
		IF_Fs = 320e3;
		rf_decim = 3;
		audio_decim = 8;
	} else if (mode == 2){
		RF_Fs = 2400e3;
		IF_Fs = 240e3;
		rf_decim = 10;
		audio_decim = 800;
		audio_expan = 147;
		BLOCK_SIZE = 1024*rf_decim*audio_decim*2;
	}
	std::vector<float> RF_h;
	std::vector<float> IF_h;
	
	std::vector<float> i_data, q_data;
	std::vector<float> filt_i, filt_q;
	std::vector<float> demod;
	std::vector<float> state_i(num_Taps, 0);
	std::vector<float> state_q(num_Taps, 0);
	std::vector<float> state_mono(num_Taps, 0);
	float prev_i = 0.0; 
	float prev_q =0.0;

//	int BLOCK_SIZE = 1024*rf_decim*audio_decim*2;
	std::vector<float> processed_data;
	auto final = 0;
	auto full_signal_start = std::chrono::high_resolution_clock::now();
	for (unsigned int block_id = 0;  ; block_id++) {
		std::vector<float> block_data(BLOCK_SIZE);
		readStdinBlockData(BLOCK_SIZE, block_id, block_data);
		if((std::cin.rdstate()) != 0) {
			std::cerr << "End of input stream reached" << std::endl;
			std::cerr << "Final run time  = "<<final<<std::endl;
			exit(1);
		}
		//std::cerr << "Read block " << block_id << std::endl;
		auto block_start = std::chrono::high_resolution_clock::now();

		split_audio_iq(block_data, i_data, q_data);

		std::cerr << "\nBlock data size: "<<block_data.size()<<std::endl;
		std::cerr << "RF Fs = "<<RF_Fs << " RF Fc = "<<RF_Fc<<std::endl;
		impulseResponseLPF(RF_Fs, RF_Fc, num_Taps, RF_h);

		//conv_ds_slow(filt_i, i_data, RF_h, rf_decim, state_i);
		//conv_ds_slow(filt_q, q_data, RF_h, rf_decim, state_q);
		// conv_ds(filt_i, i_data, RF_h, rf_decim, state_i);
		// conv_ds(filt_q, q_data, RF_h, rf_decim, state_q);
		conv_ds_fast(filt_i, i_data, RF_h, rf_decim, state_i);
		conv_ds_fast(filt_q, q_data, RF_h, rf_decim, state_q);
		FM_demod(filt_i, filt_q, prev_i, prev_q, demod);

		std::cerr << "I data size: "<< i_data.size() << std::endl;
		std::cerr << "RF H size: "<< RF_h.size()<<std::endl;
		std::cerr<< "Filtered I data size: "<< filt_i.size()<<std::endl;
		std::cerr <<"Demodulated data size: "<<demod.size()<<std::endl;
		
		if(mode == 0 || mode == 1){
			impulseResponseLPF(IF_Fs, mono_Fc, num_Taps, IF_h);

		std::cerr << "IF_h size: "<< IF_h.size() << std::endl;
		std::cerr << "IF_Fs: " << IF_Fs << " mono_Fc: "<< mono_Fc<<std::endl;
		//conv_ds_slow(processed_data, demod, IF_h, audio_decim, state_mono);
		// conv_ds(processed_data, demod, IF_h, audio_decim, state_mono);
		conv_ds_fast(processed_data, demod, IF_h, audio_decim, state_mono);
		
		} else {
			gainimpulseResponseLPF(IF_Fs*audio_expan, mono_Fc, num_Taps*audio_expan, IF_h, audio_expan);
			conv_rs(processed_data, demod, IF_h, audio_decim, audio_expan, state_mono);
		}

		std::cerr << "Read block " << block_id << " Processed_data size: " << processed_data.size() << std::endl;
		auto block_end = std::chrono::high_resolution_clock::now();
		std::chrono::duration<double, std::milli> block_time = block_end - block_start;
		std::cerr << "Block: "<< block_id<< " has runtime: "<<block_time.count()<<std::endl;
		final += block_time.count();

		std::vector<short int> audio_data(processed_data.size());
		for (unsigned int k = 0; k < processed_data.size(); k++){
			if (std::isnan(processed_data[k])) audio_data[k] = 0;
			else audio_data[k] = static_cast<short int> (processed_data[k]*16384);
		}
		fwrite(&audio_data[0], sizeof(short int),audio_data.size(),stdout);




	

	}
	auto full_signal_end = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double, std::milli> full_signal_time = full_signal_end - full_signal_start;
	std::cerr<<"Full signal took: "<<full_signal_time.count();

	return 0;
}
