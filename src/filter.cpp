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
