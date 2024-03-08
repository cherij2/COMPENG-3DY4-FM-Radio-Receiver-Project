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
    for (int i = 0, j = 0; i < input_signal.size(); i += decimation_factor, j++)
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
    for (int i = 0, j = 0; i < input_signal.size(); i += expansion_factor, j++)
    {
        output_signal[j] = input_signal[i];
    }
}
void conv_h(std::vector<float> &y, const std::vector<float> &x, const std::vector<float> &h){
    y.clear(); y.resize(x.size()+h.size());
    for(int n = 0; n<y.size(); n++){
        for(int k=0;k<(int)h.size(); k++){
            if (((n-k)>= 0) && (n-k)<(int)x.size())
            y[n] += h[k]*x[n-k];
        }
    }
}
void conv_ds_slow(std::vector<float> &y, const std::vector<float> &x, const std::vector<float> &h, int ds){
    conv_h(y, x, h);
    for(int n=0; n<int(y.size()/ds);n++){
        y[n] = y[n*ds];
    }
    y.resize(int(y.size()/ds));
    y.shrink_to_fit();
    }
    
void split_audio_iq(const std::vector<float> &audio_data, std::vector<float> &I, std::vector<float> &Q)
{
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
    for (int i = 0; i < I.size(); i++){

        float deriv_I = I[i] - I_prev;
        float deriv_Q = Q[i] - Q_prev;
        
        current_phase[i] = (I[i] == 0.0 || Q[i] == 0.0)? 0.0:(I[i] * deriv_Q - Q[i] * deriv_I) / (std::pow(I[i], 2) + std::pow(Q[i], 2));
        I_prev = I[i];
        Q_prev = Q[i];   
    }//if I or Q at that index = 0 make that demod element at that index  = 0
}
