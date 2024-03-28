#include "dy4.h"
#include "thread.h"
#include "mode.h"
#include "filter.h"
#include "RFfront.h"
#include "thread.h"
#include "rds.h"

#include <cmath>
#include <vector>

void impulseResponseRootRaisedCosine(std::vector<float> &h, float Fs, int N_taps) {
    float T_symbol = 1/2375.0;
    float beta = 0.90;
    h.clear(); h.resize(N_taps, 0.0);
    //std::vector<float> impulseResponseRRC(N_taps);

    for(int k = 0; k < N_taps; ++k) {
        float t = (k-N_taps/2)/Fs;
        if(t == 0.0) {
            h[k] = 1.0 + beta * ((4/M_PI)-1);
        } else if (t == -T_symbol / (4 * beta) || t == T_symbol / (4*beta)) {
            h[k] = (beta / std::sqrt(2)) * (((1 + 2/M_PI) * (std::sin(M_PI / (4*beta)))) + ((1 - 2/M_PI) * (std::cos(M_PI / (4*beta)))));
        } else {
            h[k] = (std::sin(M_PI * t * (1-beta) / T_symbol) +
                4*beta * t/T_symbol * std::cos(M_PI * t * (1+beta) / T_symbol)) / 
                (M_PI * t * (1 - (4 * beta * t/T_symbol) * (4*beta* t/T_symbol)) / T_symbol);
        }
    }    
}

void get_bits (const std::vector<float> &rrc_after, int sps, std::vector<float> &bits) {

    int check_within = 3;
    int offset = 0;
    int i;
    int count = 0;

    //to get first offset, change the start to the peak
    for(int k = 1; k < sps; k++) {
        // if max or min
        if((rrc_after[k] > rrc_after[k + 1] && rrc_after[k] > rrc_after[k - 1]) 
        || (rrc_after[k] < rrc_after[k + 1] && rrc_after[k] < rrc_after[k - 1])) {
            offset = k;
            count++;
        }
    }
    bits.resize(rrc_after.size(), 0);
    //this keeps running again and again and again and again and again and again 
    for(int i = offset+5; i < (int)rrc_after.size() - 5; i += sps) {
        int check_fb = i - check_within;
        //std::cerr<<"check _rb"<<std::endl;
        int check_fe = i + check_within;
        
        for(int j = check_fb; j < check_fe; j++) {
            //std::cerr<<"check _rb: "<< i << "\n" <<std::endl;
            if((rrc_after[j] > rrc_after[j + 1] && rrc_after[j] > rrc_after[j - 1]) 
            || (rrc_after[j] < rrc_after[j + 1] && rrc_after[j] < rrc_after[j - 1])){
                i = j;
            }
        }

        if(rrc_after[i] > 0) {
            //bits.push_back(1);  //local max
            bits[i] = 1;
        }  
        if(rrc_after[i] < 0) {
            //bits.push_back(-1);  //local max
            bits[i] = -1;
        }
        if(count >= std::floor(rrc_after.size() / sps)) {
            for(int t = (rrc_after.size() - sps); t < sps; t++) {
                // if max or min
                if(rrc_after[t] > 0) {
                    //bits.push_back(1);  //local max
                    bits[i] = 1;
                }  
                if(rrc_after[t] < 0) {
                    //bits.push_back(0);  //local max
                    bits[i] = -1;
                }
            }
        }
        count++;
    }
}


    
// void preDataProcessing(int mode) {
//     Mode values;
//     values.configMode(mode);

//     float RDSFb = 54000;
// 	float RDSFe = 60000;

// 	std::vector<float> rds_BPF_coeffs;
// 	std::vector<float> rds_filtered;
// 	std::vector<float> state_rds_bb(values.num_Taps-1, 0.0);

//     std::vector<float> CR_rds_BPF_coeffs;
// 	std::vector<float> CR_rds_filtered;
// 	std::vector<float> CR_state_rds(values.num_Taps-1, 0.0);

// 	std::vector<float> rds_nonlin;

// 	float CR_RDSFb = 113500;
// 	float CR_RDSFe = 11450;
//     //for BLOCK DELAY RDS
//     std::vector<float> rds_processed_delay;
//     std::vector<float> rds_state_delay((values.num_Taps-1)/2, 0.0);

//     State rds_state = {0.0, 0.0, 1.0, 0.0, 0, 1.0};
//     float rds_lockInFreq = 114000;
//     std::vector<float> rds_NCO_outp;
//     float rds_normBandwidth = 0.003;
//     float rds_phaseAdjust = 0.0;
//     float rds_ncoScale = 0.5;


//     std::vector<float> dem_mixer;
//     std::vector<float> dem_rds_LPF_coeffs;
//     std::vector<float> dem_rds_LPF_filtered;
//     std::vector<float> dem_state_rds_LPF(values.num_Taps-1, 0.0);


//     float dem_resamplerFs = 2375 * values.SPS;
//     float dem_resamplerFc = std::min((values.audio_expan / values.audio_decim) *(2375.0/2), (2375.0/2));
// 	std::vector<float> dem_rds_resamp_coeffs(values.num_Taps * values.audio_expan, 0.0);
//     std::vector<float> dem_rds_resamp_filtered;
//     std::vector<float> dem_state_rds_resamp(values.num_Taps-1, 0.0);


//     //BPF COEFFICIENTS FOR RDS BASEBAND

//     // split_audio_iq(block_data, i_data, q_data);
//     // conv_ds_fast(filt_i, i_data, RF_h, values.rf_decim, state_i);
//     // conv_ds_fast(filt_q, q_data, RF_h, values.rf_decim, state_q);
//     // FM_demod(filt_i, filt_q, prev_i, prev_q, demod);
// 	// BPFCoeffs(pilotFb, pilotFe, values.IF_Fs, values.num_Taps, pilot_BPF_coeffs);
// 	// BPFCoeffs(stereoFb, stereoFe, values.IF_Fs, values.num_Taps, stereo_BPF_coeffs);
    

//     BPFCoeffs(RDSFb, RDSFe, values.IF_Fs, values.num_Taps, rds_BPF_coeffs);
//     conv_ds_fast(rds_filtered, demod, rds_BPF_coeffs, 1, state_rds_bb);

//     // ========================
//     // CARRIER RECOVERY STARTS
//     // ========================
//     //Sqaure NonLin
//     for(int i = 0; i < rds_filtered.size(); i++) {
// 		rds_nonlin[i] = rds_filtered[i] * rds_filtered[i];
// 	} 

//     //ALL PASS FILTER
//     delayBlock(rds_filtered, rds_processed_delay, rds_state_delay);

//     //BPF CARRIER RECOVERY
//     BPFCoeffs(CR_RDSFb, CR_RDSFe, values.IF_Fs, values.num_Taps, CR_rds_BPF_coeffs);
//     conv_ds_fast(CR_rds_filtered, rds_nonlin, CR_rds_BPF_coeffs, 1, CR_state_rds);

//     //PLL CARRIER RECOVERY
//     fmPll(CR_rds_filtered, rds_NCO_outp, rds_lockInFreq, values.IF_Fs, rds_ncoScale, rds_phaseAdjust, rds_normBandwidth, rds_state);

//     // =====================
//     // CARRIER RECOVERY ENDS
//     // =====================

//     // ===================
//     // DEMODULATION STARTS
//     // ===================
    
//     //MIXER
//     for(int i = 0; i < CR_rds_filtered.size(); i++) {
// 		dem_mixer[i] += CR_rds_filtered[i] * rds_filtered[i];
// 	}

//     //LPF 
//     impulseResponseLPF(values.IF_Fs, 3000, values.num_Taps, dem_rds_LPF_coeffs);
//     conv_ds_fast(dem_rds_LPF_filtered, dem_mixer, dem_rds_LPF_coeffs, 1, dem_state_rds_LPF);

//     //CHECK IF CORRECT W/ TA, conceptually dont really get
//     //Rational Resampler ?? how does this work
    
//     impulseResponseLPF(dem_resamplerFs, dem_resamplerFc, values.num_Taps, dem_rds_resamp_coeffs);
//     //conv_rs(std::vector<float> &yb, const std::vector<float> &xb, const std::vector<float> &h, int ds, int us, std::vector<float> &state)
//     conv_rs(dem_rds_resamp_filtered, dem_rds_LPF_filtered, dem_rds_resamp_coeffs, values.audio_decim, values.audio_expan, dem_state_rds_LPF);
    
//    // impulseResponseRootRaisedCosine(float Fs, int N_taps)
// }