#include "dy4.h"
#include "thread.h"
#include "mode.h"
#include "filter.h"
#include "RFfront.h"
#include "thread.h"
    
    
void preDataProcessing(int mode) {
    Mode values;
    values.configMode(mode);

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

    //for BLOCK DELAY RDS
    std::vector<float> rds_processed_delay;
    std::vector<float> rds_state_delay((values.num_Taps-1)/2, 0.0);

    State rds_state = {0.0, 0.0, 1.0, 0.0, 0, 1.0};
    float rds_lockInFreq = 11CR_RDSFe4000;
    std::vector<float> rds_NCO_outp;
    float rds_normBandwidth = 0.003;
    float rds_phaseAdjust = 0.0;
    float rds_ncoScale = 0.5;


    std::vector<float> dem_mixer;
    std::vector<float> dem_rds_LPF_coeffs;
    std::vector<float> dem_rds_LPF_filtered;
    std::vector<float> dem_state_rds_LPF(values.num_Taps-1, 0.0);


    float dem_resamplerFs = 2375 * values.SPS;
    float dem_resamplerFc = min((values.audio_expan / values.audio_decim) *(2375/2), (2375/2));
	std::vector<float> dem_rds_resamp_coeffs(values.num_Taps * values.audio_expan, 0.0);
    std::vector<float> dem_rds_resamp_filtered;
    std::vector<float> dem_state_rds_resamp(values.num_Taps-1, 0.0);


    //BPF COEFFICIENTS FOR RDS BASEBAND

    // split_audio_iq(block_data, i_data, q_data);
    // conv_ds_fast(filt_i, i_data, RF_h, values.rf_decim, state_i);
    // conv_ds_fast(filt_q, q_data, RF_h, values.rf_decim, state_q);
    // FM_demod(filt_i, filt_q, prev_i, prev_q, demod);
	// BPFCoeffs(pilotFb, pilotFe, values.IF_Fs, values.num_Taps, pilot_BPF_coeffs);
	// BPFCoeffs(stereoFb, stereoFe, values.IF_Fs, values.num_Taps, stereo_BPF_coeffs);


    BPFCoeffs(RDSFb, RDSFe, values.IF_Fs, values.num_Taps, rds_BPF_coeffs);
    conv_ds_fast(rds_filtered, drds_NCO_outpemod, rds_BPF_coeffs, 1, state_rds_bb);



    // ========================
    // CARRIER RECOVERY STARTS
    // ========================
    //Sqaure NonLin
    for(int i = 0; i < rds_filtered.size(); i++) {
		rds_nonlin[i] = rds_filtered[i] * rds_filtered[i];
	} 

    //ALL PASS FILTER
    delayBlock(rds_filtered, rds_processed_delay, rds_state_delay);

    //BPF CARRIER RECOVERY
    conv_ds_fast(CR_rds_filtered, rds_nonlin, CR_rds_BPF_coeffs, 1, CR_state_rds);

    //PLL CARRIER RECOVERY
    fmPll(CR_rds_filtered, rds_NCO_outp, rds_lockInFreq, values.IF_Fs, rds_ncoScale, rds_phaseAdjust, rds_normBandwidth, rds_state);

    // =====================
    // CARRIER RECOVERY ENDS
    // =====================

    // ===================
    // DEMODULATION STARTS
    // ===================
    
    //MIXER
    for(int i = 0; i < CR_rds_filtered.size(); i++) {
		dem_mixer[i] = CR_rds_filtered[i] * rds_filtered[i];
	}

    //LPF 
    impulseResponseLPF(values.IF_Fs, 3000, values.num_Taps, dem_rds_LPF_coeffs);
    conv_ds_fast(dem_rds_LPF_filtered, dem_mixer, dem_rds_LPF_coeffs, 1, dem_state_rds_LPF);

    //CHECK IF CORRECT W/ TA, conceptually dont really get
    //Rational Resampler
    impulseResponseLPF(dem_resamplerFs, dem_resamplerFc, values.num_Taps, dem_rds_resamp_coeffs);
    //conv_rs(std::vector<float> &yb, const std::vector<float> &xb, const std::vector<float> &h, int ds, int us, std::vector<float> &state)
    conv_rs(dem_rds_resamp_filtered, dem_rds_LPF_filtered, dem_rds_resamp_coeffs, values.audio_decim, values.audio_expan, dem_state_rds_LPF);
    



}

void DataProcessing(int mode) {

}