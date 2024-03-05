#ifndef DY4_RFFRONT_H
#define DY4_RFFRONT_H

#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <math.h>

void downsample(const std::vector<int> &input_signal, std::vector<int> &output_signal, int decimation_factor) {}

void split_audio_iq(const std::vector<float> &audio_data, std::vector<float> &I, std::vector<float> &Q){}

void FM_demod(const std::vector<float> &I, const std::vector<float> &Q, float &prev_phase, float &I_prev, float &Q_prev, std::vector<float> &current_phase){}

#endif