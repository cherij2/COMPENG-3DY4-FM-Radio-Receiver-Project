

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

	int mode = 0;
	int channel = 0;
	if(argc < 2){
		std::cerr<< "Operating in default mode 0" << std::endl;
	}
	else if (argc == 2){

		mode = atoi(argv[1]);
		if (mode > 3){
			std::cerr << "Wrong mode " << mode << std::endl;
			exit(1);		
		}
		else {
			std::cerr << "Usage: "<< argv[0] <<std::endl;
			std::cerr << "or " << std::endl;
			std::cerr << "Usage: "<< argv[0] << "<mode>" << std::endl;
			std::cerr << "\t\t <mode> is a value from 0 to 3" << std::endl;
			exit(1);
			}
		std::cerr << "Operating in mode" << mode << std::endl;
	}
	else if (argc == 3){
		mode = atoi(argv[1]);
		channel = atoi(argv[2]);
		if (mode > 3){
			std::cerr << "Wrong mode "<< mode << std::endl;
			exit(1);
		}
		else {
			std::cerr << "Usage: "<< argv[0] <<std::endl;
			std::cerr << "or " << std::endl;
			std::cerr << "Usage: "<< argv[0] << "<mode>" << std::endl;
			std::cerr << "\t\t <mode> is a value from 0 to 3" << std::endl;
			exit(1);
			}
		}
	} */
	
	int main(int argc, char* argv[])
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
