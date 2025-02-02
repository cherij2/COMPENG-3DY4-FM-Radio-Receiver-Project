#ifndef DY4_MODE_H
#define DY4_MODE_H

struct Mode{
    int RF_Fs;
	float RF_Fc;
	float IF_Fs;
	float mono_Fc;
	float num_Taps; //if too high of a value takes too long to run
	unsigned short int rf_decim;
	float audio_decim;
	float audio_expan;
	int BLOCK_SIZE;
	int SPS;
	int rds_up;
	int rds_down;

Mode();

void configMode(int mode);


};

#endif