#include "dy4.h"
#include "RFfront.h"

void downsample(const std::vector<int> &input_signal, std::vector<int> &output_signal, int decimation_factor)
{
    // Calculate the size of the downsampled vector
    int downsampled_size = (input_signal.size() + decimation_factor - 1) / decimation_factor;

    // Resize the output signal vector to hold the downsampled values
    output_signal.clear();
    output_signal.resize(downsampled_size);

    // Perform the downsampling
    for (int i = 0, j = 0; i < (int)input_signal.size(); i += decimation_factor, j++)
    {
        output_signal[j] = input_signal[i];
    }
}
void upsample(const std::vector<int> &input_signal, std::vector<int> &output_signal, int expansion_factor)
{
    // Calculate the size of the downsampled vector
    int upsampled_size = (input_signal.size() + expansion_factor - 1) / expansion_factor;

    // Resize the output signal vector to hold the downsampled values
    output_signal.clear();
    output_signal.resize(upsampled_size);

    // Perform the downsampling
    for (int i = 0, j = 0; i < (int)input_signal.size(); i += expansion_factor, j++)
    {
        output_signal[j] = input_signal[i];
    }
}

void blockConvolutionFIR(std::vector<float> &yb, const std::vector<float> &xb, const std::vector<float> &h, std::vector<float> &state)
{				// parameters include yb which is output block, xb input signal, h is the impulse response of the filter, state is state that will be used and updated to be used for the next block
	yb.clear(); // this implementation copies the python code
	yb.resize(xb.size(), 0.0);

	for (int n = 0; n < (int)xb.size(); n++)
	{
		for (int k = 0; k < (int)h.size(); k++)
		{
			if (n - k >= 0)
			{
				yb[n] += h[k] * xb[n - k];
			}
			else
			{
				yb[n] += h[k] * state[h.size() - 1 + (n - k)]; //
			}
		}
	}

	for (int i = 0; i < (int)h.size() - 1; i++) // updates the state at the end for the next block to be used
	{
		state[i] = xb[xb.size() - h.size() + 1 + i];
	}
}

void conv_h(std::vector<float>& y, const std::vector<float>& x, const std::vector<float>& h) {
    y.clear(); y.resize(x.size() + h.size());
    for (int n = 0; n < (int)y.size(); n++) {
        for (int k = 0; k < (int)h.size(); k++) {
            if (((n - k) >= 0) && ((n - k) < (int)x.size())) {
                y[n] += h[k] * x[n - k];
            }
        }
    }
}

void conv_ds_slow(std::vector<float>& y, const std::vector<float>& x, const std::vector<float>& h, int ds, std::vector<float> &state) {
    conv_h(y, x, h);
    //blockConvolutionFIR(y, x, h, state);
    for (int n = 0; n < int(y.size() / ds); n++) {
        y[n] = y[n * ds];
    }
    y.resize(int(y.size() / ds));
    y.shrink_to_fit();
}

void conv_ds(std::vector<float> &yb, const std::vector<float> &xb, const std::vector<float> &h, int ds, std::vector<float> &state){				// parameters include yb which is output block, xb input signal, h is the impulse response of the filter, state is state that will be used and updated to be used for the next block
	yb.clear(); // this implementation copies the python code
	yb.resize(xb.size()/ds, 0.0);

	for (int n = 0; n < (int)xb.size(); n++)
	{
        float sum = 0.0;
		for (int k = 0; k < (int)h.size(); k++)
		{
			if (n - k >= 0)
			{
				sum += h[k] * xb[n - k];
			}
			else
			{
				sum += h[k] * state[state.size() - 1 + (n - k)]; //
			}

            if(n%ds == 0){
                if(n == 0){
                    yb[0] = sum;
                }
                else{
                    yb[n/ds] = sum;
                }
            }
		}
	}
    std::vector<float> new_state(&xb[xb.size()-state.size()], &xb[xb.size()]);
    state = new_state;
    //std::cerr<<"new_state of "<<new_state.size()<<std::endl;
}
void conv_ds_fast(std::vector<float> &yb, const std::vector<float> &xb, const std::vector<float> &h, int ds, std::vector<float> &state){				// parameters include yb which is output block, xb input signal, h is the impulse response of the filter, state is state that will be used and updated to be used for the next block
	yb.clear(); // this implementation copies the python code
	yb.resize(xb.size()/ds, 0.0);

	for (int n = 0; n < (int)yb.size(); n++)
	{
		for (int k = 0; k < (int)h.size(); k++)
		{
			if (ds*n - k >= 0)
			{
				yb[n] += h[k] * xb[ds*n - k];
			}
			else
			{
				yb[n] += h[k] * state[state.size() - 1 + (ds*n - k)]; //
			}

		}
	}
    std::vector<float> new_state(&xb[xb.size()-state.size()], &xb[xb.size()]);
    state = new_state;
    //std::cerr<<"new_state: "<<new_state.size()<<std::endl;
}


void conv_rs(std::vector<float> &yb, const std::vector<float> &xb, const std::vector<float> &h, int ds, int us, std::vector<float> &state){
    yb.clear(); yb.resize((xb.size()*us)/ds, 0.0);
    //fast implementation from lecture notes
    //std::cerr<< "yb: "<<yb.size()<<" xb size: "<<xb.size()<<" h size: "<<h.size()<<std::endl;
    
    for(int n = 0; n < (int)yb.size(); n++){
        int phase = (n*ds)%us; // phase changes when n is incremented in our case
        for (int k = phase; k < (int)h.size(); k+= us){
            int dx = (ds*n-k)/us; //when n = 1, dx starts at 2, 1,  
            
            if(dx >= 0){
                yb[n] += h[k]*xb[dx];
                
            }else{
                 yb[n] += h[k]*state[state.size() - 1 +dx];
                 
            }
        }
    }
    //std::cerr<< "resampling state size = "<<state.size()<<std::endl;
    std::vector<float> new_state(&xb[xb.size()-state.size()], &xb[xb.size()]);
    //std::cerr<<"new_state size: "<<new_state.size()<<std::endl;
    state = new_state;
}

void split_audio_iq(const std::vector<float> &audio_data, std::vector<float> &I, std::vector<float> &Q)
{
    I.clear(); //I.resize(audio_data.size()/2);
    Q.clear(); //Q.resize(audio_data.size(/2);
    for (int i = 0; i < (int)audio_data.size(); i++)
    {
        if (i % 2 == 0)
            I.push_back(audio_data[i]);
        else
            Q.push_back(audio_data[i]);
    }
}

void FM_demod(const std::vector<float> &I, const std::vector<float> &Q, float &I_prev, float &Q_prev, std::vector<float> &current_phase)
{

    current_phase.clear(); current_phase.resize(I.size());
    for (int i = 0; i < (int)I.size(); i++){

        float denominator = (std::pow(I[i], 2) + std::pow(Q[i], 2));
        float deriv_I = I[i] - I_prev;
        float deriv_Q = Q[i] - Q_prev;
        // if (I[i] == 0 || Q[i] == 0 || denominator == 0){
        //     current_phase[i] = 0;
        // }
        // else {
        current_phase[i] = (I[i] == 0.0 || Q[i] == 0.0) ? 0.0:(I[i] * deriv_Q - Q[i] * deriv_I) / denominator;
        //above just implements an if else statement and assigns the value of demod at index i 0 if 

        // } 
        //std::cerr<<" FM demod at index: "<<i<<" FM demod value: "<<current_phase[i]<<std::endl;
        I_prev = I[i];
        Q_prev = Q[i];   
    }//if I or Q at that index = 0 make that demod element at that index  = 0
}
