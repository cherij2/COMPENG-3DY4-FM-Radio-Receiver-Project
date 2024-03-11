#ifndef DY4_RFFRONT_H
#define DY4_RFFRONT_H

#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <math.h>

void downsample(const std::vector<int> &input_signal, std::vector<int> &output_signal, int decimation_factor);

void split_audio_iq(const std::vector<float> &audio_data, std::vector<float> &I, std::vector<float> &Q);

void FM_demod(const std::vector<float> &I, const std::vector<float> &Q, float &I_prev, float &Q_prev, std::vector<float> &current_phase);

void blockConvolutionFIR(std::vector<float> &yb, const std::vector<float> &xb, const std::vector<float> &h, std::vector<float> &state);

void conv_h(std::vector<float> &y, const std::vector<float> &x, const std::vector<float> &h);

void conv_ds_slow(std::vector<float> &y, const std::vector<float> &x, const std::vector<float> &h, int ds, std::vector<float> &state);

void conv_ds(std::vector<float> &yb, const std::vector<float> &xb, const std::vector<float> &h, int ds, std::vector<float> &state);

#endif