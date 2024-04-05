#include "mode.h"


Mode::Mode(): RF_Fs(2400e3), RF_Fc(100e3), IF_Fs(240e3), mono_Fc(16e3), num_Taps(101),
        rf_decim(10), audio_decim(5), audio_expan(1), BLOCK_SIZE(1024*rf_decim*audio_decim*2), SPS (30), rds_up(19.0), rds_down(64.0) {}

void Mode::configMode(int mode){
    switch (mode){
        case 0:
        //default is set in the constructor
            break;
        case 1:
            RF_Fs = 960e3;
            IF_Fs = 320e3;
            rf_decim = 3;
            audio_decim = 8;
            audio_expan = 1;
            BLOCK_SIZE = 1500*rf_decim*audio_decim;
            SPS = 1;
            rds_up = 1.0;
            rds_down = 1.0;
            break;
        case 2:
            RF_Fs = 2400e3;
            IF_Fs = 240e3;
            rf_decim = 10;
            audio_decim = 800;
            audio_expan = 147;
            BLOCK_SIZE = audio_decim*audio_expan;
            SPS = 43;
            rds_up = 817.0;
            rds_down = 1920.0;
            break;
        case 3:
            RF_Fs = 960e3;
            IF_Fs = 120e3;
            rf_decim = 8;
            audio_decim = 400;
            audio_expan = 147;
            BLOCK_SIZE = 12*audio_decim*rf_decim*2;
            RF_Fc = 60e3;
            SPS = 1;
            rds_up = 1.0;
            rds_down = 1.0;
            break;
        }
    }
