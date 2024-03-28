#include "mode.h"


Mode::Mode(): RF_Fs(2400e3), RF_Fc(100e3), IF_Fs(240e3), mono_Fc(16e3), num_Taps(101),
        rf_decim(10), audio_decim(5), audio_expan(1), BLOCK_SIZE(1024*rf_decim*audio_decim*2){}

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
            BLOCK_SIZE = 800*rf_decim*audio_decim*2;
            break;
        case 2:
            RF_Fs = 2400e3;
            IF_Fs = 240e3;
            rf_decim = 10;
            audio_decim = 800;
            audio_expan = 147;
            BLOCK_SIZE = 5*audio_decim*rf_decim*2;
            break;
        case 3:
            RF_Fs = 960e3;
            IF_Fs = 120e3;
            rf_decim = 8;
            audio_decim = 400;
            audio_expan = 147;
            BLOCK_SIZE = 15*audio_decim*rf_decim*2;
            RF_Fc = 60e3;
            break;
        }
    }
