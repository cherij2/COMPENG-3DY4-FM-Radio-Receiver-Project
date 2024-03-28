

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
#include "mode.h"
#include "thread.h"
#include "rds.h"


int main(int argc, char *argv[])
{

	// Default mode 0
	int mode = 0;
	std::string channel = "m";

	// Mode Selection
	if (argc<2){
		std::cerr << "Operating in default mode 0" << std::endl;
	} else if (argc==2){
		mode=atoi(argv[1]);
		if (mode>3){
			std::cerr << "Wrong mode " << mode << std::endl;
			exit(1);
		}
	} else if(argc == 3){
		mode = atoi(argv[1]);
		channel = argv[2];
		if (mode > 3){
			std::cerr << "Wrong mode "<< mode << std::endl;
		} else if (channel != "m" || channel != "s" || channel != "r"){
			std::cerr << "wrong channel "<<channel << std::endl;
		}
	} 
	else {
		std::cerr << "Usage: " << argv[0] << std::endl;
		std::cerr << "or " << std::endl;
		std::cerr << "Usage: " << argv[0] << " <mode>" << std::endl;
		std::cerr << "\t\t <mode> is a value from 0-3" << argv[0] << std::endl;
		exit(1);
	}
	std::cerr << "Operating in mode " << mode << std::endl;

	// for (int i = 0; i < 233; i++ ){
	// 	std::cerr<<"index i: "<<i<<std::endl;
	// 	std::vector<float>final_output = produce_data(mode);
	// }
	// rf_thread(mode);
	// audio_thread(mode);
	std::thread rf_producer(rf_thread, mode);  // Create the RF producer thread
    std::thread audio_consumer(audio_thread, mode);  // Create the audio consumer thread

    // Wait for both threads to finish
    
	rf_producer.join();
    audio_consumer.join();
	
	//--------------------------------UNCOMMENT BELOW TO RUN RDS WITHOUT THREADING-------------------
	/*
	Mode values;
	values.configMode(mode);
	
	
	int RF_Fs = values.RF_Fs;
	float RF_Fc = values.RF_Fc;
	float IF_Fs = values.IF_Fs;
	float mono_Fc = values.mono_Fc;
	float num_Taps = values.num_Taps; //if too high of a value takes too long to run
	unsigned short int rf_decim = values.rf_decim;
	float audio_decim = values.audio_decim;
	float audio_expan = values.audio_expan;
	int BLOCK_SIZE = values.BLOCK_SIZE;
	float RDSFb = 54000;
	float RDSFe = 60000;

	std::vector<float> rds_BPF_coeffs;
	std::vector<float> rds_filtered;
	std::vector<float> state_rds_bb(values.num_Taps-1, 0.0);

    std::vector<float> CR_rds_BPF_coeffs;
	std::vector<float> CR_rds_filtered;
	std::vector<float> CR_state_rds(values.num_Taps-1, 0.0);

	std::vector<float> rds_nonlin;

	float CR_RDSFb = 113500;
	float CR_RDSFe = 114500;

	//WHAT DO FOR MODE 1 AND 3??
	float rds_fs_rat_res = values.SPS * 2375;
	std::vector<float> outp_rrc;
	std::vector<float> state_rrc(values.num_Taps-1, 0.0);


    //for BLOCK DELAY RDS
    std::vector<float> rds_processed_delay;
    std::vector<float> rds_state_delay((values.num_Taps-1)/2, 0.0);

    State rds_state = {0.0, 0.0, 1.0, 0.0, 0, 1.0, 1.0};
    float rds_lockInFreq = 114000;
    std::vector<float> rds_NCO_outp;
	std::vector<float> rds_NCO_outpQ;
    float rds_normBandwidth = 0.003;
    float rds_phaseAdjust = 0.0;
    float rds_ncoScale = 0.5;


    std::vector<float> dem_mixer;
    std::vector<float> dem_rds_LPF_coeffs;
    std::vector<float> dem_rds_LPF_filtered;
    std::vector<float> dem_state_rds_LPF(values.num_Taps-1, 0.0);


    float dem_resamplerFs = 2375 * values.SPS;
    float dem_resamplerFc = std::min((values.audio_expan / values.audio_decim) *(2375.0/2), (2375.0/2));
	std::vector<float> dem_rds_resamp_coeffs(values.num_Taps * values.audio_expan, 0.0);
    std::vector<float> dem_rds_resamp_filtered;
    std::vector<float> dem_state_rds_resamp(values.num_Taps-1, 0.0);

	std::vector<float> RRC_coeffs;
	std::vector<float> RF_h;
	// std::vector<float> final_coeffs;

	std::vector<float> i_data, q_data;
	std::vector<float> filt_i, filt_q;
	std::vector<float> demod;
	std::vector<float> state_i(num_Taps-1, 0.0);
	std::vector<float> state_q(num_Taps-1, 0.0);
	// std::vector<float> state_mono(num_Taps-1, 0.0);

	// float pilotFb = 18500;
	// float pilotFe = 19500;
	// float stereoFb = 22000;
	// float stereoFe = 54000;

	// float STnumTaps = 101;



	// std::vector<float> pilot_BPF_coeffs;
	// std::vector<float> stereo_BPF_coeffs;
	// std::vector<float> pilot_filtered;
	// std::vector<float> stereo_filtered;
	// std::vector<float> state_pilot(num_Taps-1, 0.0);
	// std::vector<float> state_stereo(num_Taps-1, 0.0);



	// //MIXER VARIABLES
	// std::vector<float> mixer;
	// std::vector<float> mixer_coeffs;
	// std::vector<float> mixer_filtered;
	// std::vector<float> state_mixer(num_Taps-1, 0.0);

	// //for BLOCK DELAY
	// std::vector<float> mono_processed_delay;
	// std::vector<float> state_delay((num_Taps-1)/2, 0.0);

	// std::vector<float> right_stereo;
	// std::vector<float> left_stereo;

	// State state = {0.0, 0.0, 1.0, 0.0, 0, 1.0};
	// float pilot_lockInFreq = 19000;
	// std::vector<float> pilot_NCO_outp;
	// float normBandwidth = 0.01;
	// float phaseAdjust = 0.0;
	// float ncoScale = 2.0;


	// std::cerr<<"TEST"<<std::endl;
	float prev_i = 0.0;
	float prev_q =0.0;
	// std::vector<float> processed_data;
	// std::vector<float> stereo_data;
	// std::cerr<<" num_taps "<<values.num_Taps<<std::endl;
	// //LPF COEFFICIENTS FOR FRONT END AND MONO PATH
	gainimpulseResponseLPF(RF_Fs, RF_Fc, num_Taps, RF_h, audio_expan); //FRONT END
	// gainimpulseResponseLPF(IF_Fs*audio_expan, mono_Fc, num_Taps*audio_expan, final_coeffs, audio_expan);//MONO PATH
	BPFCoeffs(RDSFb, RDSFe, values.IF_Fs, values.num_Taps, rds_BPF_coeffs);
	//BPF CARRIER RECOVERY
    BPFCoeffs(CR_RDSFb, CR_RDSFe, values.IF_Fs, values.num_Taps, CR_rds_BPF_coeffs);
	// //BPF COEFFICIENTS FOR STEREO PILOT FREQUENCY 1ST and STEREOBAND 2ND
	// BPFCoeffs(pilotFb, pilotFe, IF_Fs, STnumTaps, pilot_BPF_coeffs);
	// BPFCoeffs(stereoFb, stereoFe, IF_Fs, STnumTaps, stereo_BPF_coeffs);
	impulseResponseLPF(values.IF_Fs, 3000, values.num_Taps, dem_rds_LPF_coeffs);
	//impulseResponseLPF(dem_resamplerFs, dem_resamplerFc, values.num_Taps, dem_rds_resamp_coeffs);
	gainimpulseResponseLPF(dem_resamplerFs, dem_resamplerFc, values.num_Taps, dem_rds_resamp_coeffs, values.rds_up);
	// float final = 0.0;//THIS HOLDS THE FINAL RUN TIME OF MONO PATH FOR NOW
	// auto full_signal_start = std::chrono::high_resolution_clock::now();
	while (true){
		for (unsigned int block_id = 0;  ; block_id++) {
			std::vector<float> block_data(BLOCK_SIZE);
			readStdinBlockData(BLOCK_SIZE, block_id, block_data); //block_data holds the data for one block
			//if((std::cin.rdstate()) != 0) {
			if(block_id == 21){//RUNS UNTIL BLOCK 21 THEN BREAKS (USED FOR TESTING)
				std::cerr << "End of input stream reached" << std::endl;
				//FINAL RUN TIME IS THE ADDITION OF THE RUN TIME FOR EACH BLOCK
				//std::cerr << "Final run time  = "<<final<<std::endl;
				exit(1);
			}

			//--------------------RF-FRONT END-----------------------
			//STD CERR WAS USED FOR DEBUGGING MOST OF THE ISSUES
			//std::cerr << "Read block " << block_id << std::endl;
			//auto block_start = std::chrono::high_resolution_clock::now();
			//std::cerr<<"Mono Cutoff: "<<mono_Fc<<std::endl;
			std::cerr<<"before split iq"<<std::endl;
			split_audio_iq(block_data, i_data, q_data);

			//std::cerr << "\nBlock data size: "<<block_data.size()<<std::endl;
			//std::cerr << "RF Fs = "<<RF_Fs << " RF Fc = "<<RF_Fc<<std::endl;
			//COULD IMPLEMENT A FUNCTION THAT DOES CONVOLUTION FOR I AND Q IN ONE RUN

			conv_ds_fast(filt_i, i_data, RF_h, rf_decim, state_i);
			std::cerr<<" test "<<std::endl;
			conv_ds_fast(filt_q, q_data, RF_h, rf_decim, state_q);
			FM_demod(filt_i, filt_q, prev_i, prev_q, demod);
			std::cerr<<"after fm demod"<<std::endl;
			

	// 		//--------------------END OF RF-FRONT END-------------------

	//		//---------------------START OF RDS-------------------------
			// std::cerr<<"before rds"<<std::endl;
			conv_ds_fast(rds_filtered, demod, rds_BPF_coeffs, 1, state_rds_bb);
			// std::cerr<<"after rds 0"<<std::endl;
			rds_nonlin.resize(rds_filtered.size());
			for(int i = 0; i < rds_filtered.size(); i++) {
				rds_nonlin[i] = rds_filtered[i] * rds_filtered[i];
			}
			// std::cerr<<"after non linearity"<<std::endl;
			//ALL PASS FILTER
    		delayBlock(rds_filtered, rds_processed_delay, rds_state_delay);
			conv_ds_fast(CR_rds_filtered, rds_nonlin, CR_rds_BPF_coeffs, 1, CR_state_rds);
			// std::cerr<<"tst"<<std::endl;
			//PLL
			fmPll(CR_rds_filtered, rds_NCO_outp, rds_NCO_outpQ, rds_lockInFreq, values.IF_Fs, rds_ncoScale, rds_phaseAdjust, rds_normBandwidth, rds_state); 		
			// std::cerr<<"after pll"<<std::endl;


			dem_mixer.resize(CR_rds_filtered.size());
			for(int i = 0; i<CR_rds_filtered.size();i++){
				dem_mixer[i] = 3*rds_NCO_outp[i] * rds_processed_delay[i];
			}
			
			conv_ds_fast(dem_rds_LPF_filtered, dem_mixer, dem_rds_LPF_coeffs, 1, dem_state_rds_LPF);

    		conv_rs(dem_rds_resamp_filtered, dem_rds_LPF_filtered, dem_rds_resamp_coeffs, values.rds_down, values.rds_up, dem_state_rds_LPF);
			
			impulseResponseRootRaisedCosine(RRC_coeffs, rds_fs_rat_res, values.num_Taps);
			conv_ds_fast(outp_rrc, dem_rds_resamp_filtered, RRC_coeffs, 1, state_rrc);
			std::vector<float> bits;
			// std::cerr<<"before get bits"<<std::endl;
			get_bits(outp_rrc, values.SPS, bits);
			// std::cerr<<"after getting bits"<<std::endl;

			
			std::cerr<<" dem resamp filtered size "<<dem_rds_resamp_filtered.size()<<std::endl;
			if (block_id == 20){
				std::vector<float> pllinput_index;
				genIndexVector(pllinput_index, CR_rds_filtered.size());
				logVector("CR_rds_filtered", pllinput_index, CR_rds_filtered);
				std::vector<float> vector_index;
				genIndexVector(vector_index, outp_rrc.size());
				logVector("outp_rrc", vector_index, outp_rrc);
				std::vector<float> dem_mixer_index;
				genIndexVector(dem_mixer_index, rds_NCO_outp.size());
				logVector("rds_NCO_outp", dem_mixer_index, rds_NCO_outp);
				std::vector<float>output_bits;
				genIndexVector(output_bits, bits.size());
				logVector("bits", output_bits, bits);
				std::vector<float> pre_cdr_q;
				genIndexVector(pre_cdr_q, rds_NCO_outpQ.size());
				logVector("rds_NCO_outpQ", pre_cdr_q, rds_NCO_outpQ);
			}
	
		}//ends for


		return 0;
	}//ends while
	*/
}
