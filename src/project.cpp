

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

	// DEFAULT VALUES ASSUMING MODE 0
	int RF_Fs = 2400e3;
	float RF_Fc = 100e3;
	float IF_Fs = 240e3;
	float mono_Fc;
	float num_Taps = 101; //if too high of a value takes too long to run
	unsigned short int rf_decim = 10;
	float audio_decim = 5;
	float audio_expan = 1;
	int BLOCK_SIZE = 1024*rf_decim*audio_decim*2;

	//MODE SELECT
	if (mode == 1){ //only downsampling in this mode
		RF_Fs = 960e3;
		IF_Fs = 320e3;
		rf_decim = 3;
		audio_decim = 8;
		audio_expan = 1;
		BLOCK_SIZE = 1500*audio_decim*rf_decim*2;
	} else if (mode == 2){ //resampling needed in this mode
		RF_Fs = 2400e3;
		IF_Fs = 240e3;
		rf_decim = 10;
		audio_decim = 800;
		audio_expan = 147;
		BLOCK_SIZE = audio_decim*audio_expan;
	} else if(mode == 3){ //resampling needed in this mode
		RF_Fs = 960e3;
		IF_Fs = 120e3;
		rf_decim = 8;
		audio_decim = 400;
		audio_expan = 147;
		BLOCK_SIZE = 15*audio_decim*rf_decim*2;
	}
	std::cerr<<"audio expan: "<<audio_expan<<" audio decim: "<<audio_decim<<std::endl;
	std::cerr<<"min between these two: " <<(audio_expan/audio_decim)<<" "<<(IF_Fs/2)<<" and this valie: "<<IF_Fs/2<<std::endl;
	mono_Fc = ((std::min((int)((audio_expan/audio_decim)*(IF_Fs/2)), (int)IF_Fs/2)) < 16000) ? (std::min((int)((audio_expan/audio_decim)*(IF_Fs/2)), (int)IF_Fs/2)) : 16000.0;

	std::vector<float> RF_h;
	std::vector<float> final_coeffs;

	std::vector<float> i_data, q_data;
	std::vector<float> filt_i, filt_q;
	std::vector<float> demod;
	std::vector<float> state_i(num_Taps, 0.0);
	std::vector<float> state_q(num_Taps, 0.0);
	std::vector<float> state_mono(num_Taps, 0.0);

	float pilotFb = 18500;
	float pilotFe = 19500;
	float stereoFb = 22000;
	float stereoFe = 54000;

	float STnumTaps = 101;

	float errorD = 0.0;
	float integrator = 0.0;
	float phaseEst = 0.0;
	float feedbackI = 1.0;
	float feedbackQ = 0.0;
	int trigOffset = 0;

	std::vector<float> pilot_BPF_coeffs;
	std::vector<float> stereo_BPF_coeffs;
	std::vector<float> pilot_filtered;
	std::vector<float> stereo_filtered;
	std::vector<float> state_pilot(num_Taps-1, 0.0);
	std::vector<float> state_stereo(num_Taps-1, 0.0);



	//MIXER VARIABLES
	std::vector<float> mixer;
	std::vector<float> mixer_coeffs;
	std::vector<float> mixer_filtered;
	std::vector<float> state_mixer(num_Taps-1, 0.0);

	//for BLOCK DELAY
	std::vector<float> mono_processed_delay;
	std::vector<float> state_delay((num_Taps-1)/2, 0.0);

	std::vector<float> right_stereo;
	std::vector<float> left_stereo;

	State state = {0.0, 0.0, 1.0, 0.0, 0, 1.0};
	float pilot_lockInFreq = 19000;
	std::vector<float> pilot_NCO_outp;
	float normBandwidth = 0.01;
	float phaseAdjust = 0.0;
	float ncoScale = 2.0;


	std::cerr<<"TEST";
	float prev_i = 0.0;
	float prev_q =0.0;
	std::vector<float> processed_data;
	std::vector<float> stereo_data;
	//LPF COEFFICIENTS FOR FRONT END AND MONO PATH
	gainimpulseResponseLPF(RF_Fs, RF_Fc, num_Taps, RF_h, audio_expan); //FRONT END
	gainimpulseResponseLPF(IF_Fs*audio_expan, mono_Fc, num_Taps*audio_expan, final_coeffs, audio_expan);//MONO PATH

	//BPF COEFFICIENTS FOR STEREO PILOT FREQUENCY 1ST and STEREOBAND 2ND
	BPFCoeffs(pilotFb, pilotFe, IF_Fs, STnumTaps, pilot_BPF_coeffs);
	BPFCoeffs(stereoFb, stereoFe, IF_Fs, STnumTaps, stereo_BPF_coeffs);
	

	float final = 0.0;//THIS HOLDS THE FINAL RUN TIME OF MONO PATH FOR NOW
	auto full_signal_start = std::chrono::high_resolution_clock::now();
	while (true){
		for (unsigned int block_id = 0;  ; block_id++) {
			std::vector<float> block_data(BLOCK_SIZE);
			readStdinBlockData(BLOCK_SIZE, block_id, block_data); //block_data holds the data for one block
			if((std::cin.rdstate()) != 0) {
			//if(block_id == 4){
				std::cerr << "End of input stream reached" << std::endl;
				//FINAL RUN TIME IS THE ADDITION OF THE RUN TIME FOR EACH BLOCK
				std::cerr << "Final run time  = "<<final<<std::endl;
				exit(1);
			}
			//--------------------RF-FRONT END-----------------------
			//STD CERR WAS USED FOR DEBUGGING MOST OF THE ISSUES
			//std::cerr << "Read block " << block_id << std::endl;
			auto block_start = std::chrono::high_resolution_clock::now();
			//std::cerr<<"Mono Cutoff: "<<mono_Fc<<std::endl;

			split_audio_iq(block_data, i_data, q_data);

			//std::cerr << "\nBlock data size: "<<block_data.size()<<std::endl;
			//std::cerr << "RF Fs = "<<RF_Fs << " RF Fc = "<<RF_Fc<<std::endl;
			//COULD IMPLEMENT A FUNCTION THAT DOES CONVOLUTION FOR I AND Q IN ONE RUN

			conv_ds(filt_i, i_data, RF_h, rf_decim, state_i);
			conv_ds(filt_q, q_data, RF_h, rf_decim, state_q);
			FM_demod(filt_i, filt_q, prev_i, prev_q, demod);

			//--------------------END OF RF-FRONT END-------------------

			// std::cerr << "I data size: "<< i_data.size() << std::endl;
			// std::cerr << "RF H size: "<< RF_h.size()<<std::endl;
			// std::cerr<< "Filtered I data size: "<< filt_i.size()<<std::endl;
			// std::cerr <<"Demodulated data size: "<<demod.size()<<std::endl;


			delayBlock(demod, mono_processed_delay, state_delay);
			//-------------------MONO PATH START------------------------
			//WE CAN USE THE RESAMPLING FUNCTION BECAUSE WE ASSIGN audio_expan a value of 1

			// std::cerr << "final_coeffs size: "<< final_coeffs.size() << std::endl;
			// std::cerr << "IF_Fs: " << IF_Fs << " mono_Fc: "<< mono_Fc<<std::endl;

			conv_rs(processed_data, mono_processed_delay, final_coeffs, audio_decim, audio_expan, state_mono);

			//-------------------MONO PATH END--------------------------

			//-------------------STEREO PATH START--------------------------
			// fmPll(const std::vector<float>& pllIn, std::vector<float>& ncoOut, float freq, float Fs, float ncoScale = 1.0, float phaseAdjust = 0.0, float normBandwidth = 0.01)
			//bandPass(float Fb, float Fe, float Fs, unsigned short int num_taps, std::vector &h)


			//to get pilot freq

			std::cerr<<"pilot BPF coeffs: "<<pilot_BPF_coeffs.size()<<std::endl;
			// for (int i = 0; i<pilot_BPF_coeffs.size(); i++){
			// 	std::cerr<<"BPF within 18.5k and 19.5k: "<<pilot_BPF_coeffs[i]<<std::endl;


			// }
			// for (int i = 0; i<pilot_BPF_coeffs.size(); i++){
			// std::cerr<<"LPF for IF: "<<final_coeffs[i]<<std::endl;}
			conv_ds_fast(pilot_filtered, demod, pilot_BPF_coeffs, 1, state_pilot);

			//convolveFIR(std::vector<float> &y, const std::vector<float> &x, const std::vector<float> &h)

			//to get stereo band freq
		
			conv_ds_fast(stereo_filtered, demod, stereo_BPF_coeffs, 1, state_stereo);
			// for (int i = 0; i<stereo_BPF_coeffs.size(); i++){
			// 	std::cerr<<"BPF within 22k and 54k: "<<stereo_BPF_coeffs[i]<<std::endl;


			// }
			std::cerr<<"pilot_BPF_coeffs size: "<< pilot_BPF_coeffs.size()<<"\tstereo_BPF_coeffs: "<<stereo_BPF_coeffs.size()<<std::endl;
			std::cerr<<"pilot filtered size: "<<pilot_filtered.size()<< "\tstereo_filtered size: "<<stereo_filtered.size()<<std::endl;




			// for (int i = 0; i < 30;i++){
			// 	std::cerr<<"pilot filtered at "<<i<<" "<<pilot_filtered[i]<<std::endl;
			// }
			// for(int i = 0; i <30;i++){
			// 	std::cerr<<"pilot band coeffs at "<<i<<" "<< pilot_BPF_coeffs[i]<<std::endl;
			// }
			//PLL for pilot
			//void fmPll(const std::vector<float>& pllIn, std::vector<float>& ncoOut, float freq, float Fs, float integrator, float phaseEst, float feedbackI, float feedbackQ, int trigOffset, float ncoScale = 2.0, float phaseAdjust = 0.0, float normBandwidth = 0.01)
	
			fmPll(pilot_filtered, pilot_NCO_outp, pilot_lockInFreq, IF_Fs, ncoScale, phaseAdjust, normBandwidth, state);
			 


			//MIXER
			mixer.resize(stereo_filtered.size(), 0.0);
			for(int i = 0; i < stereo_filtered.size(); i++) {
				mixer[i] = 2 * pilot_NCO_outp[i] * stereo_filtered[i];
				// if(i < 30){
				// std::cerr<< "index: "<<i<<"\tpilot NCO output "<<pilot_NCO_outp[i]<<"\tstereo filtered"<<stereo_filtered[i]<<"\tmixer val: "<<mixer[i]<<std::endl;
				
				// }
			}

			

			conv_ds_fast(mixer_filtered, mixer, final_coeffs, audio_decim, state_mixer);
			// for (int i = 0; i < 30; i++){
			// 	std::cerr<<"mixer coeffs: "<<final_coeffs[i]<<std::endl;
			// }

			//RIGHT AND LEFT STEREO
			right_stereo.resize(mixer_filtered.size());
			left_stereo.resize(mixer_filtered.size());
			for(int i = 0; i < mixer_filtered.size(); i++) {
				//!!!! is equation correct?
				right_stereo[i] = (mixer_filtered[i] - processed_data[i]);
				left_stereo[i] = (mixer_filtered[i] + processed_data[i]);
				// if(i < 30){
				// std::cerr<<"Mixer filtered at index: "<<i<<" value: "<<mixer_filtered[i]<<std::endl;
				// }
			}



			std::cerr<<"Mixer filtered size: "<<mixer_filtered.size()<<std::endl;


			//figure out how to implement 'state saving' in
			//finish implementing delay function in filter
			stereo_data.resize(right_stereo.size()*2);
			int i = 0;
			for (int k = 0; k< right_stereo.size(); k++){
				stereo_data[i] = left_stereo[k];
				stereo_data[i+1] = right_stereo[k];
				i += 2;
			
			}
			std::cerr<<"Stereo Data size: "<<stereo_data.size()<<"Left and right: "<<left_stereo.size()<<right_stereo.size()<<std::endl;
			for (int i  = 0; i < 21;i++){
				std::cerr<<"Stereo data at index: "<<i<<" data: "<<stereo_data[i]<<std::endl;
				std::cerr<<"Right data: "<<right_stereo[i]<<std::endl;
				std::cerr<<"Left data: "<<left_stereo[i]<<std::endl;
			}

			//-------------------STEREO PATH END--------------------------


			std::cerr << "Read block " << block_id << " Processed_data size: " << processed_data.size() << std::endl;

			//BELOW SHOWS THE RUN TIME FOR EACH BLOCK AFTER CONVOLUTION IS RUN
			auto block_end = std::chrono::high_resolution_clock::now();
			std::chrono::duration<double, std::milli> block_time = block_end - block_start;
			std::cerr << "Block: "<< block_id<< " has runtime: "<<block_time.count()<<"\n"<<std::endl;
			final += block_time.count();
			///------------BELOW WRITES THE MONO PATH------------
			//CODE BELOW IS WHAT WRITES THE AUDIO IF NAN assigns audio at k = 0;
			// std::vector<short int> audio_data(processed_data.size());
			// for (unsigned int k = 0; k < processed_data.size(); k++){
			// 	if (std::isnan(processed_data[k])) audio_data[k] = 0;
			// 	else audio_data[k] = static_cast<short int> (processed_data[k]*16384); //MULTIPLYING BY 16384 NORMALIZES DATA B/W -1 and 1

			// }
			// //WRITES AUDIO TO STANDARD OUTPUT AS 16 bit
			// fwrite(&audio_data[0], sizeof(short int),audio_data.size(),stdout);

			//------------BELOW WRITES ONLY ONE CHANNEL THE LEFT------------
			// std::vector<short int> audio_data(left_stereo.size());
			// for (unsigned int k = 0; k < left_stereo.size(); k++){
			// 	if (std::isnan(left_stereo[k])) audio_data[k] = 0;
			// 	else audio_data[k] = static_cast<short int> (left_stereo[k]*16384); //MULTIPLYING BY 16384 NORMALIZES DATA B/W -1 and 1
			// }
			// //WRITES AUDIO TO STANDARD OUTPUT AS 16 bit
			// fwrite(&audio_data[0], sizeof(short int),audio_data.size(),stdout);

			//------------BELOW WRITES ONLY ONE CHANNEL THE RIGHT------------
			std::vector<short int> audio_data(right_stereo.size());
			for (unsigned int k = 0; k < right_stereo.size(); k++){
				if (std::isnan(right_stereo[k])) audio_data[k] = 0;
				else audio_data[k] = static_cast<short int> (right_stereo[k]*16384); //MULTIPLYING BY 16384 NORMALIZES DATA B/W -1 and 1
			}
			//WRITES AUDIO TO STANDARD OUTPUT AS 16 bit
			fwrite(&audio_data[0], sizeof(short int),audio_data.size(),stdout);

		}//ends for


		return 0;
	}//ends while
}
