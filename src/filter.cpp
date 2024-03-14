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


void fmPll(const std::vector<float>& pllIn, std::vector<float>& ncoOut, float &freq, float &Fs, float &integrator, float &phaseEst, float &feedbackI, float &feedbackQ, int &trigOffset, float &errorD, float ncoScale = 2.0, float phaseAdjust = 0.0, float normBandwidth = 0.01) {
    // Constants for the loop filter
    const float Cp = 2.666;
    const float Ci = 3.555;

    // Calculate gains
    const float Kp = normBandwidth * Cp;
    const float Ki = normBandwidth * normBandwidth * Ci;

    // Resize and initialize output vector
    ncoOut.resize(pllIn.size() + 1);
    ncoOut[0] = 1.0;

    // Loop through input samples
    for (size_t k = 0; k < pllIn.size(); ++k) {
        // Phase detector
        float errorI = pllIn[k] * feedbackI; // In-phase error
        float errorQ = pllIn[k] * (-feedbackQ); // Quadrature error

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
    }
}

