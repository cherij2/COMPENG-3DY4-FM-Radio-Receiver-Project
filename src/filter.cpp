/*
Comp Eng 3DY4 (Computer Systems Integration Project)

Department of Electrical and Computer Engineering
McMaster University
Ontario, Canada
*/

#include "dy4.h"
#include "filter.h"

// function to compute the impulse response "h" based on the sinc function
void impulseResponseLPF(float Fs, float Fc, unsigned short int num_taps, std::vector<float> &h){
	// allocate memory for the impulse response
	h.clear();
	h.resize(num_taps, 0.0);
	float Norm_cutoff = Fc / (Fs / 2);

	for (int m = 0; m < num_taps; m++)
	{
		if (m == ((num_taps - 1) / 2))
		{
			h[m] = Norm_cutoff;
		}
		else
		{
			h[m] = Norm_cutoff * (std::sin(PI * Norm_cutoff * (m - (num_taps - 1) / 2)) / (PI * Norm_cutoff * (m - (num_taps - 1) / 2)));
		}
		h[m] = h[m] * pow((std::sin((m * PI) / num_taps)), 2);
	}
}

void gainimpulseResponseLPF(float Fs, float Fc, unsigned short int num_taps, std::vector<float> &h, int U){
	// allocate memory for the impulse response
	h.clear();
	h.resize(num_taps, 0.0);
	float Norm_cutoff = Fc / (Fs / 2);

	for (int m = 0; m < num_taps; m++)
	{
		if (m == ((num_taps - 1) / 2))
		{
			h[m] = Norm_cutoff;
		}
		else
		{
			h[m] = Norm_cutoff * (std::sin(PI * Norm_cutoff * (m - (num_taps - 1) / 2)) / (PI * Norm_cutoff * (m - (num_taps - 1) / 2)));
		}
		h[m] = U*h[m] * pow((std::sin((m * PI) / num_taps)), 2);
	}
}
// function to compute the filtered output "y" by doing the convolution
// of the input data "x" with the impulse response "h"
void convolveFIR(std::vector<float> &y, const std::vector<float> &x, const std::vector<float> &h){
	// allocate memory for the output (filtered) data
	y.clear();
	y.resize(x.size() + h.size() - 1, 0.0);
	for (int n = 0; n < (int)y.size(); n++)
	{
		y[n] = 0;
		for (int k = 0; k < (int)h.size(); k++)
		{
			if (n - k >= 0 && n - k < (int)x.size())
			{
				y[n] += h[k] * x[n - k];
			}
		}
	}
	// the rest of the code in this function is to be completed by yo	// based on your understanding and the Python code from the first lab
}


void BPFCoeffs(float Fb, float Fe, float Fs, unsigned short int num_taps, std::vector<float> &h) {
	h.clear();
	h.resize(num_taps, 0.0);
	float norm_center = ((Fe+Fb)/2)/(Fs/2) ;
	float norm_pass = (Fe-Fb)/(Fs/2) ; 
	
	for(int i = 0; i < num_taps - 1; i++) {
		if(i == ((num_taps-1)/2)) {
			h[i] = norm_pass;
		} else {
			h[i] = norm_pass * ((sin(PI*(norm_pass/2)*(i-(num_taps-1)/2)))/(PI*(norm_pass/2)*(i-(num_taps-1)/2)));
		}
		
		h[i] = h[i]*cos((i-(num_taps-1)/2)*PI*norm_center);
		h[i] = h[i]*sin(i*PI/num_taps)*sin(i*PI/num_taps);
	}
}


void fmPll(const std::vector<float>& pllIn, std::vector<float>& ncoOut, float freq, float Fs, float &integrator, float &phaseEst, float &feedbackI, float &feedbackQ, int &trigOffset, float &errorD, float ncoScale, float phaseAdjust, float normBandwidth) {
    // Constants for the loop filter
    const float Cp = 2.666;
    const float Ci = 3.555;
	std::cerr<< "norm Bandwidth: "<<normBandwidth<<"\tphaseAdjust: "<<phaseAdjust<<"\tncoScale: "<<ncoScale<<std::endl;
    // Calculate gains
    const float Kp = normBandwidth * Cp;
    const float Ki = normBandwidth * normBandwidth * Ci;

    // Resize and initialize output vector
    ncoOut.resize(pllIn.size() + 1);
    ncoOut[0] = ;

    // Loop through input samples
    for (size_t k = 0; k < pllIn.size(); ++k) {
        // Phase detector
        float errorI = pllIn[k] * feedbackI; // Infreqor
        float errorQ = pllIn[k] * (-feedbackQ); //freqe error

        // Arc tangent phase discriminator
        float errorD = std::atan2(errorQ, errorI);

        // Loop filter
        integrator += Ki * errorD;

        // Update phase estimate
        phaseEst += Kp * errorD + integrator;

        // Update internal oscillator state
        trigOffset++;

        float trigArg = 2 * M_PI * (freq / Fs) * trigOffset + phaseEst;
        feedbackI = std::cos(trigArg);
        feedbackQ = std::sin(trigArg);
        ncoOut[k + 1] = std::cos(trigArg * ncoScale + phaseAdjust);
		//std::cerr<<"ncdOut at index"<< k << " + 1 "<<ncoOut[k+1]<<std::endl;
    }
	std::cerr<<"ncoOUT at [0]: "<<ncoOut[0]<<"\tncoOut[5220]"<<ncoOut[pllIn.size()]<<std::endl;
	
	std::cerr<<"ncoOUT size: "<<ncoOut.size()<<std::endl;
}

void delayBlock(std::vector<float> input_block, std::vector<float> &output_block, int num_taps, std::vector<float> &state) {
	output_block.clear();
	output_block.resize(input_block.size());
	state.resize(num_taps/2);
	
	output_block.assign(state.begin(), state.end());
	output_block.insert(output_block.end(), input_block.begin(), input_block.end() - state.size());
	state.assign(input_block.end() - state.size(), input_block.end());
}

// void blockConvolutionFIR(std::vector<float> &yb, const std::vector<float> &xb, const std::vector<float> &h, std::vector<float> &state)
// {				// parameters include yb which is output block, xb input signal, h is the impulse response of the filter, state is state that will be used and updated to be used for the next block
// 	yb.clear(); // this implementation copies the python code
// 	yb.resize(xb.size(), 0.0);

// 	for (int n = 0; n < (int)xb.size(); n++)
// 	{
// 		for (int k = 0; k < (int)h.size(); k++)
// 		{
// 			if (n - k >= 0)
// 			{
// 				yb[n] += h[k] * xb[n - k];
// 			}
// 			else
// 			{
// 				yb[n] += h[k] * state[h.size() - 1 + (n - k)]; //
// 			}
// 		}
// 	}

// 	for (int i = 0; i < (int)h.size() - 1; i++) // updates the state at the end for the next block to be used
// 	{
// 		state[i] = xb[xb.size() - h.size() + 1 + i];
// 	}
// }