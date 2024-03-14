/*
Comp Eng 3DY4 (Computer Systems Integration Project)

Department of Electrical and Computer Engineering
McMaster University
Ontario, Canada
*/

#ifndef DY4_FILTER_H
#define DY4_FILTER_H

// add headers as needed
#include <iostream>
#include <fstream>
#include <vector>
#include <iomanip>
#include <complex>
#include <cmath>
#include <stdlib.h>
#include <algorithm>

// declaration of a function prototypes
void impulseResponseLPF(float, float, unsigned short int, std::vector<float> &);
void gainimpulseResponseLPF(float Fs, float Fc, unsigned short int num_taps, std::vector<float> &h, int U);
void convolveFIR(std::vector<float> &, const std::vector<float> &, const std::vector<float> &);
void blockConvolutionFIR(std::vector<float> &yb, const std::vector<float> &xb, const std::vector<float> &h, std::vector<float> &state);
void BPFCoeffs(float Fb, float Fe, float Fs, unsigned short int num_taps, std::vector<float> &h);
void fmPll(const std::vector<float>& pllIn, std::vector<float>& ncoOut, float freq, float Fs, float &integrator, float &phaseEst, float &feedbackI, float &feedbackQ, int &trigOffset, float &errorD, float ncoScale = 2.0, float phaseAdjust = 0.0, float normBandwidth = 0.01);

#endif // DY4_FILTER_H
