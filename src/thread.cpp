#include "dy4.h"
#include "thread.h"
#include "mode.h"
#include "filter.h"
#include "RFfront.h"

std::atomic<bool> done{false};
std::atomic<float> total_block_time;


template <typename T>
class ThreadSafeQueue {
private:
    mutable std::mutex m;                    // Mutex to protect access to the queue
    std::queue<std::shared_ptr<T>> q;        // Standard queue wrapped inside the thread-safe queue
    std::condition_variable cv;              // Condition variable for notifying waiting threads
    std::condition_variable cv_producer;     // Condition variable for notifying when space is available
    const size_t max_size = 4;               // Maximum size of the queue

public:
    // Enqueues an element by copying it into a shared_ptr and adding it to the queue
    void push(T value) {
        std::shared_ptr<T> data(std::make_shared<T>(std::move(value)));
        std::unique_lock<std::mutex> lock(m); // Use unique_lock to be able to wait

        // Wait until there is space in the queue
        cv_producer.wait(lock, [this]{ return q.size() < max_size; });

        q.push(data);                        // Push the data onto the queue
        cv.notify_one();                     // Notify one waiting thread that there is new data available
    }

    // Waits for an element to be available and pops it from the queue
    std::shared_ptr<T> wait_and_pop() {
        std::unique_lock<std::mutex> lock(m); // Unique lock allows the lock to be temporarily released
        cv.wait(lock, [this]{ return !q.empty(); }); // Wait until the queue is not empty
        std::shared_ptr<T> result = q.front(); // Get the front element
        q.pop();                              // Remove the element from the queue
        cv_producer.notify_one();             // Notify one waiting producer that space is available
        return result;                        // Return the data to the caller
    }    

    // Checks if the queue is empty
    bool empty() const {
        //std::lock_guard<std::mutex> lock(m);  // Lock the mutex to synchronize access to the empty check
        return q.empty();                     // Return whether the queue is empty
    }

    size_t size() const {
        //std::lock_guard<std::mutex> lock(m);
        return q.size(); // Get the size of the queue
    }

    // Method to print the addresses of the contents of the queue
    void print_contents() const {
        std::lock_guard<std::mutex> lock(m);
        std::queue<std::shared_ptr<T>> temp_queue = q; // Make a copy of the queue

        std::cout << "Queue contents' addresses: [";
        while (!temp_queue.empty()) {
            const auto& item_ptr = temp_queue.front(); // Access the shared_ptr to the front element
            temp_queue.pop(); // Remove the element from the temporary queue

            if (item_ptr) { // Check if the shared_ptr actually points to an object
                // Print the address pointed to by the shared_ptr
                std::cout << static_cast<const void*>(item_ptr.get()) << " ";
            } else {
                std::cout << "nullptr "; // In case the shared_ptr is null
            }

            if (!temp_queue.empty()) {
                std::cout << ", "; // Print a comma unless it's the last element
            }
        }
        std::cout << "]" << std::endl;
    }
};


// ==================================

ThreadSafeQueue<std::vector<float>> tsQueue; // Global instance of the thread-safe queuedemod_ptr

// // Function representing the work of the RF thread (the producer)
void rf_thread(int mode)  {                        // Continue producing until done is true
    Mode values;
    values.configMode(mode);
    std::vector<float> i_data, q_data;
	std::vector<float> filt_i, filt_q;
    std::vector<float> state_i(values.num_Taps-1, 0.0);
	std::vector<float> state_q(values.num_Taps-1, 0.0);
    std::vector<float> RF_h;
    std::vector<float> demod;
    float prev_i = 0.0;
    float prev_q = 0.0;
    float final_block_time = 0.0;
    float final_split_time = 0.0;
    float final_conv_i_time = 0.0;
    float final_conv_q_time = 0.0;
    float final_demod_time = 0.0;
    float final_enqueue_time = 0.0;
    float final_read_block_in = 0.0;
    auto impulse_start = std::chrono::high_resolution_clock::now();
    gainimpulseResponseLPF(values.RF_Fs, values.RF_Fc, values.num_Taps, RF_h, values.audio_expan);
    auto impulse_end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> impulse_time = impulse_end - impulse_start;
    float impulse_time_block = impulse_time.count();
    std::cerr << "NUMBER OF TAPS = "<<values.num_Taps<<std::endl;
    bool exitwhile = false;
    while (!exitwhile){
        for(unsigned int block_id = 0; ; block_id++){
            total_block_time = 0.0;
            auto block_start = std::chrono::high_resolution_clock::now();
            std::cerr<<"Reading block id "<<block_id<<std::endl;
            std::vector<float> block_data(values.BLOCK_SIZE);
            auto read_block_start = std::chrono::high_resolution_clock::now();
            readStdinBlockData(values.BLOCK_SIZE, block_id, block_data);
            auto read_block_end = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double, std::milli> read_block_time = read_block_end - read_block_start;
            final_read_block_in += read_block_time.count();
            if((std::cin.rdstate()) != 0){
            //if(block_id == 100){
                std::cerr << "BLOCK SIZE "<<values.BLOCK_SIZE<<std::endl;
                std::cerr<<"End of input stream reached" << std::endl;
                std::cerr << "FRONT END IMPULSE TIME = "<<impulse_time_block<< " ms"<<std::endl;
                std::cerr << "NUMBER OF TAPS = "<<values.num_Taps << std::endl;
                std::cerr << "RUNTIME OF SPLITTING IQ = "<<final_split_time << " ms"<<std::endl;
                std::cerr << "RUNTIME OF CONV I = "<<final_conv_i_time << " ms"<<std::endl;
                std::cerr << "RUNTIME OF CONV Q = "<<final_conv_q_time << " ms"<<std::endl;
                std::cerr << "RUNTIME OF FINAL DEMOD = "<<final_demod_time<< " ms"<<std::endl;
                std::cerr << "RUNTIME OF ENQUEUE OPERATION = "<<final_enqueue_time<< " ms"<<std::endl;
                std::cerr << "RUNTIME OF WHOLE BLOCK = "<<final_block_time<< " ms"<<std::endl;
                //tsQueue.print_contents();
                std::cerr<<"size of queue: "<<tsQueue.size()<<std::endl;
                exitwhile = true;
                std::cerr<<"done flag "<<done<<std::endl;
                done = true;
                std::cerr<<"done flag "<<done<<std::endl;
                break;
            }
            //auto function_start = std::chrono::high_resolution_clock::now();
            // std::cerr<<"Block id "<<block_id<<std::endl;
            auto split_start = std::chrono::high_resolution_clock::now();
            split_audio_iq(block_data, i_data, q_data);
            auto split_end = std::chrono::high_resolution_clock::now();
            auto conv_i_start = std::chrono::high_resolution_clock::now();
            conv_ds_fast(filt_i, i_data, RF_h, values.rf_decim, state_i);
            auto conv_i_end = std::chrono::high_resolution_clock::now();
            auto conv_q_start = std::chrono::high_resolution_clock::now();
            conv_ds_fast(filt_q, q_data, RF_h, values.rf_decim, state_q);
            auto conv_q_end = std::chrono::high_resolution_clock::now();
            auto demod_start = std::chrono::high_resolution_clock::now();
            FM_demod(filt_i, filt_q, prev_i, prev_q, demod);
            auto demod_end = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double, std::milli> split_time = split_end - split_start;
            std::chrono::duration<double, std::milli> conv_i_time = conv_i_end - conv_i_start;
            std::chrono::duration<double, std::milli> conv_q_time = conv_q_end - conv_q_start;
            std::chrono::duration<double, std::milli> demod_time = demod_end - demod_start;
            final_split_time += split_time.count();
            final_conv_i_time += conv_i_time.count();
            final_conv_q_time += conv_q_time.count();
            final_demod_time += demod_time.count();
            //---------------------UNCOMMENT FOR RDS THREADING 
            //auto demod_data = std::make_shared<std::vector<float>>(demod);
            //std::cerr<<"size before pushing: "<<tsQueue.size()<<std::endl;
            //tsQueue.push(demod_data);
            auto queue_push_start = std::chrono::high_resolution_clock::now();
            tsQueue.push(demod);
            auto queue_push_end = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double, std::milli> queue_push_time = queue_push_end - queue_push_start;
            final_enqueue_time += queue_push_time.count();
            //auto block_end = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double, std::milli> block_time = queue_push_end - block_start;
            final_block_time += block_time.count();
            std::cerr << "Block: "<< block_id << " has runtime: " << block_time.count() << std::endl;
            //total_block_time  = total_block_time + block_time.count();
            //std::cerr<<"size after pushing: "<<tsQueue.size()<<std::endl;
            std::cerr<<"\n";
        }
    }
        // Push the produced data onto the queue

}


// // Function representing the work of the audio thread (the consumer)
void audio_thread(int mode, std::string channel) {
    Mode values;
    values.configMode(mode);
    int block_count = 0;
    //std::cerr<<"ENTERED AUDIO THREAD"<<std::endl;
    std::vector<float> state_mono(values.num_Taps-1, 0.0);
    std::vector<float> state_stereo(values.num_Taps-1, 0.0); // band pass 22k-54k
    std::vector<float> state_pilot(values.num_Taps-1, 0.0); // band pass 18.5khz-19.5khz (noise one)
    std::vector<float> state_mixer(values.num_Taps-1, 0.0); // multiple stereo bands passes floats together

    std::vector<float> pilot_BPF_coeffs;
	std::vector<float> stereo_BPF_coeffs;
	std::vector<float> pilot_filtered;
	std::vector<float> stereo_filtered;

	// //MIXER VARIABLES
	std::vector<float> mixer;
	std::vector<float> mixer_coeffs;
	std::vector<float> mixer_filtered;

	// //for BLOCK DELAY
	std::vector<float> mono_processed_delay;
	std::vector<float> state_delay((values.num_Taps-1)/2, 0.0);

	std::vector<float> right_stereo;
	std::vector<float> left_stereo;

	State state = {0.0, 0.0, 1.0, 0.0, 0, 1.0};
	float pilot_lockInFreq = 19000;
	std::vector<float> pilot_NCO_outp;
    std::vector<float> pilot_NCO_outpQ;
	float normBandwidth = 0.01;
	float phaseAdjust = 0.0;
	float ncoScale = 2.0;

    float pilotFb = 18500;
	float pilotFe = 19500;
	float stereoFb = 22000;
	float stereoFe = 54000;

	std::vector<float> processed_data;
	std::vector<float> stereo_data;

    std::vector<float> final_coeffs;

    bool exitwhile = false;

    float final_mono_impulse = 0.0;
    auto mono_impulse_start = std::chrono::high_resolution_clock::now();
	gainimpulseResponseLPF(values.IF_Fs*values.audio_expan, values.mono_Fc, values.num_Taps*values.audio_expan, final_coeffs, values.audio_expan);//MONO PATH
    auto mono_impulse_end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> mono_impulse_time = mono_impulse_end - mono_impulse_start;
    final_mono_impulse += mono_impulse_time.count();
    auto bpf_pilot_start = std::chrono::high_resolution_clock::now();
	BPFCoeffs(pilotFb, pilotFe, values.IF_Fs, values.num_Taps, pilot_BPF_coeffs);
    auto bpf_pilot_end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> bpf_pilot_time = bpf_pilot_end - bpf_pilot_start;
    auto bpf_extract_start = std::chrono::high_resolution_clock::now();
	BPFCoeffs(stereoFb, stereoFe, values.IF_Fs, values.num_Taps, stereo_BPF_coeffs);
    auto bpf_extract_end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> bpf_extract_time = bpf_extract_end - bpf_extract_start;
    
    float final_pop_time = 0.0;
    float final_mono_convrs_time = 0.0;
    float final_mono_path_time = 0.0;
    //float final_delay_time = 0.0;
    float final_stereo_channel_time = 0.0;
    float final_stereo_carrier_recovery_time = 0.0;
    float final_stereo_processing_time = 0.0;
    float final_stereo_path_time = 0.0;
    float final_delay_time = 0.0;
    float final_conv_pilot_time = 0.0;
    float final_conv_stereo_channel_time = 0.0;
    float final_pll_time = 0.0;
    float final_mixer_time = 0.0;
    float final_stereo_conv_rs_time = 0.0;
    float final_lr_channel_time = 0.0;
    float final_interleave_time = 0.0;

    
    while (!exitwhile) {                           // Continue consuming until done is true
        auto mono_block_start = std::chrono::high_resolution_clock::now();
        auto stereo_block_start = std::chrono::high_resolution_clock::now();
        if(tsQueue.empty() && done) { // if queue is empty and there is nothing else that is coming in, then break
            std::cerr << "FINAL MONO IMPULSE RESPONSE GENERATION = "<<final_mono_impulse<< " ms"<<std::endl;
            std::cerr << "FINAL POP TIME = "<<final_pop_time<< " ms"<<std::endl;
            std::cerr << "FINAL MONO RESAMPLING = "<<final_mono_convrs_time<< " ms"<<std::endl;
            std::cerr << "FINAL MONO RUNTIME = "<<final_mono_path_time<< " ms"<<std::endl;
            std::cerr << "FINAL STEREO CHANNEL EXTRACTION TIME = "<<final_stereo_channel_time<< " ms"<<std::endl;
            std::cerr << "FINAL STEREO CARRIER RECOVERY TIME = "<<final_stereo_carrier_recovery_time<< " ms"<<std::endl;
            std::cerr << "FINAL STEREO PROCESSING TIME = "<<final_stereo_processing_time<< " ms"<<std::endl;
            std::cerr << "FINAL STEREO PATH TIME = "<<final_stereo_path_time<< " ms"<<std::endl;
            std::cerr << "FINAL DELAY TIME = "<<final_delay_time<< " ms"<<std::endl;
            std::cerr << "FINAL CONV PILOT = "<<final_conv_pilot_time<< " ms"<<std::endl;
            std::cerr << "FINAL CONV STEREO TIME = "<<final_conv_stereo_channel_time<< " ms"<<std::endl;
            std::cerr << "FINAL PLL TIME = "<<final_pll_time<< " ms"<<std::endl;
            std::cerr << "FINAL MIXER TIME = "<<final_mixer_time<< " ms"<<std::endl;
            std::cerr << "FINAL STEREO RESAMPLING = "<<final_stereo_conv_rs_time<< " ms"<<std::endl;
            std::cerr << "FINAL LR CHANNEL = "<<final_lr_channel_time<< " ms"<<std::endl;
            std::cerr << "FINAL FINAL INTERLEAVE TIME = "<<final_interleave_time<<" ms"<<std::endl;
            exitwhile = true;
            break;
        }
        //int i  = 0;
        //std::cerr<<"test"<<std::endl;
        auto pop_start = std::chrono::high_resolution_clock::now();
        std::shared_ptr<std::vector<float>> demod_ptr = tsQueue.wait_and_pop(); // Wait for and pop data from the queue
        auto pop_end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double, std::milli> pop_time = pop_end - pop_start;
        final_pop_time += pop_time.count();
        //std::cerr << "Processing block " << block_count << std::endl;
        // std::cerr<<"i val: "<<i<<"demod ptr "<<demod_ptr<<std::endl;
        // i++;
        // if (demod_ptr) {
        //     // Dereference the pointer to obtain the vector
        //     const std::vector<float>& demod_vector = *demod_ptr;

        //     // Print out the contents of the vector
        //     std::cerr<<"demid vector size "<<demod_vector.size()<<std::endl;
        // }
        //-------------MONO PATH START--------------------------------
        auto delay_start = std::chrono::high_resolution_clock::now();
        delayBlock(*demod_ptr, mono_processed_delay, state_delay);
        auto delay_end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double, std::milli> delay_time = delay_end - delay_start;
        final_delay_time += delay_time.count();
        //std::cerr<<"mono process delay size "<<mono_processed_delay.size()<<std::endl;
        auto mono_convrs_start = std::chrono::high_resolution_clock::now();
        conv_rs(processed_data, mono_processed_delay, final_coeffs, values.audio_decim, values.audio_expan, state_mono);
        auto mono_convrs_end = std::chrono::high_resolution_clock::now();
        auto mono_block_end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double, std::milli> mono_convrs_time = mono_convrs_end - mono_convrs_start;
        std::chrono::duration<double, std::milli> mono_path_time = mono_block_end - mono_block_start;
        final_mono_path_time += mono_path_time.count();
        final_mono_convrs_time += mono_convrs_time.count();
        final_delay_time += delay_time.count();

        //std::cerr<<"processed data size: "<<processed_data.size()<<std::endl;
        //-------------MONO PATH END----------------------------------

        //fmPll(const std::vector<float>& pllIn, std::vector<float>& ncoOut, float freq, float Fs, float ncoScale = 1.0, float phaseAdjust = 0.0, float normBandwidth = 0.01)
        auto stereo_channel_start = std::chrono::high_resolution_clock::now();
        conv_ds_fast(stereo_filtered, *demod_ptr, stereo_BPF_coeffs, 1, state_stereo);
        auto stereo_channel_end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double, std::milli> stereo_channel_time = stereo_channel_end - stereo_channel_start;
        final_stereo_channel_time += stereo_channel_time.count()+bpf_extract_time.count();
        final_conv_stereo_channel_time += stereo_channel_time.count();

        auto stereo_carrier_start = std::chrono::high_resolution_clock::now();
        conv_ds_fast(pilot_filtered, *demod_ptr, pilot_BPF_coeffs, 1, state_pilot);
        auto stereo_conv_pilot_end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double, std::milli> conv_pilot_time = stereo_conv_pilot_end - stereo_carrier_start;
        final_conv_pilot_time += conv_pilot_time.count();
        // std::cerr<<"pilot filtered size: "<<pilot_filtered.size()<<"stereo filtered size: "<<stereo_filtered.size()<<std::endl;
        auto fm_pll_start = std::chrono::high_resolution_clock::now();
        fmPll(pilot_filtered, pilot_NCO_outp,pilot_NCO_outpQ, pilot_lockInFreq, values.IF_Fs, ncoScale, phaseAdjust, normBandwidth, state);
        auto fm_pll_end = std::chrono::high_resolution_clock::now();
        auto stereo_carrier_end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double, std::milli> pll_time = fm_pll_end - fm_pll_start;
        final_pll_time += pll_time.count();
        std::chrono::duration<double, std::milli> stereo_carrier_recovery_time = stereo_carrier_end - stereo_carrier_start;
        final_stereo_carrier_recovery_time += stereo_carrier_recovery_time.count()+bpf_pilot_time.count();
        auto stereo_processing_start = std::chrono::high_resolution_clock::now();
        auto mixer_start = std::chrono::high_resolution_clock::now();
		mixer.resize(stereo_filtered.size(), 0.0);
		for(unsigned int i = 0; i < stereo_filtered.size(); i++) {
            mixer[i] = 2 * pilot_NCO_outp[i] * stereo_filtered[i];
        }
        auto mixer_end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double, std::milli> mixer_time = mixer_end - mixer_start;
        final_mixer_time += mixer_time.count();


        auto stereo_conv_rs_start = std::chrono::high_resolution_clock::now();
        conv_rs(mixer_filtered, mixer, final_coeffs, values.audio_decim, values.audio_expan, state_mixer);
        auto stereo_conv_rs_end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double, std::milli> stereo_conv_rs_time = stereo_conv_rs_end - stereo_conv_rs_start;
        final_stereo_conv_rs_time += stereo_conv_rs_time.count();
        // std::cerr<<"mixer size: "<<mixer.size()<<" mixer filtered size: "<<mixer_filtered.size()<<std::endl;
        auto lr_channel_start = std::chrono::high_resolution_clock::now();
        right_stereo.resize(mixer_filtered.size());
        left_stereo.resize(mixer_filtered.size());
        for(unsigned int i = 0; i < mixer_filtered.size(); i++) {
            right_stereo[i] = (mixer_filtered[i] - processed_data[i]);
            left_stereo[i] = (mixer_filtered[i] + processed_data[i]);
        }
        auto lr_channel_end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double, std::milli> lr_channel_time = lr_channel_end - lr_channel_start;
        final_lr_channel_time += lr_channel_time.count();

        // std::cerr<<"left size: "<<left_stereo.size()<<" right stereo size: "<<right_stereo.size()<<std::endl;
        // std::cerr<<"right stereo size" <<right_stereo.size()<<std::endl;
        // std::cerr<<"left stereo size" <<left_stereo.size()<<std::endl;
        auto interleave_start = std::chrono::high_resolution_clock::now();
        stereo_data.resize(right_stereo.size()*2);
        int i = 0;
        for (unsigned int k = 0; k< right_stereo.size(); k++){
            stereo_data[i] = left_stereo[k];
            stereo_data[i+1] = right_stereo[k];
            i += 2;
        }
        auto interleave_end = std::chrono::high_resolution_clock::now();
        auto stereo_processing_end = std::chrono::high_resolution_clock::now();
        auto stereo_block_end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double, std::milli> stereo_processing_time = stereo_processing_end - stereo_processing_start;
        final_stereo_processing_time += stereo_processing_time.count();
        std::chrono::duration<double, std::milli> stereo_path_time = stereo_block_end - stereo_block_start;
        final_stereo_path_time += stereo_path_time.count();
        std::chrono::duration<double, std::milli> interleave_time = interleave_end - interleave_start;
        final_interleave_time += interleave_time.count();

            //WRITE MONO CHANNEL AUDIO TO STANDARD OUTPUT AS 16 BIT
            if (channel == "m"){
            std::vector<short int> audio_data(processed_data.size());
                for (unsigned int k = 0; k < processed_data.size(); k++){
                    if (std::isnan(processed_data[k])) audio_data[k] = 0;
                    else audio_data[k] = static_cast<short int> (processed_data[k]*16384); //MULTIPLYING BY 16384 NORMALIZES DATA B/W -1 and 1

                }
                //WRITES STEREO CHANNEL AUDIO TO STANDARD OUTPUT AS 16 bit
                fwrite(&audio_data[0], sizeof(short int),audio_data.size(),stdout);
            } else if(channel == "s"){
                std::vector<short int> audio_data(stereo_data.size());
                for (unsigned int k = 0; k < stereo_data.size(); k++){
                    if (std::isnan(stereo_data[k])) audio_data[k] = 0;
                    else audio_data[k] = static_cast<short int> (stereo_data[k]*16384); //MULTIPLYING BY 16384 NORMALIZES DATA B/W -1 and 1
                }
                //WRITES RIGHT CHANNEL AUDIO TO STANDARD OUTPUT AS 16 bit
                fwrite(&audio_data[0], sizeof(short int),audio_data.size(),stdout);
            } else if (channel == "r"){
                std::vector<short int> audio_data(right_stereo.size());
			    for (unsigned int k = 0; k < right_stereo.size(); k++){
				    if (std::isnan(right_stereo[k])) audio_data[k] = 0;
				    else audio_data[k] = static_cast<short int> (right_stereo[k]*16384); //MULTIPLYING BY 16384 NORMALIZES DATA B/W -1 and 1
			}
			//WRITES LEFT CHANNEL AUDIO TO STANDARD OUTPUT AS 16 bit
			fwrite(&audio_data[0], sizeof(short int),audio_data.size(),stdout);
            } else if (channel == "l"){
                std::vector<short int> audio_data(left_stereo.size());
			    for (unsigned int k = 0; k < left_stereo.size(); k++){
				    if (std::isnan(left_stereo[k])) audio_data[k] = 0;
				    else audio_data[k] = static_cast<short int> (left_stereo[k]*16384); //MULTIPLYING BY 16384 NORMALIZES DATA B/W -1 and 1
			}
			//WRITES AUDIO TO STANDARD OUTPUT AS 16 bit
			fwrite(&audio_data[0], sizeof(short int),audio_data.size(),stdout);
            }
        //total_block_time = .count();
        //std::cerr << "BLOCK RUNTIME: "<<total_block_time<<std::endl;
    //
        block_count++;
    }
}
//-----------------------RDS THREAD BELOW---------------------
// void rdsThread(int mode) {
// 	//--------------RDS INITIALIZATION------------------
//Mode values;
//    values.configMode(mode);
// 	float RDSFb = 54000;
// 	float RDSFe = 60000;

// 	std::vector<float> rds_BPF_coeffs;
// 	std::vector<float> rds_filtered;
// 	std::vector<float> state_rds_bb(values.num_Taps-1, 0.0);

//     std::vector<float> CR_rds_BPF_coeffs;
// 	std::vector<float> CR_rds_filtered;
// 	std::vector<float> CR_state_rds(values.num_Taps-1, 0.0);

// 	std::vector<float> rds_nonlin;

// 	float CR_RDSFb = 113500;
// 	float CR_RDSFe = 114500;

// 	//WHAT DO FOR MODE 1 AND 3??
// 	float rds_fs_rat_res = values.SPS * 2375;
// 	std::vector<float> outp_rrc;
// 	std::vector<float> state_rrc(values.num_Taps-1, 0.0);


//     //for BLOCK DELAY RDS
//     std::vector<float> rds_processed_delay;
//     std::vector<float> rds_state_delay((values.num_Taps-1)/2, 0.0);

//     State rds_state = {0.0, 0.0, 1.0, 0.0, 0, 1.0, 1.0};
//     float rds_lockInFreq = 114000;
//     std::vector<float> rds_NCO_outp;
// 	std::vector<float> rds_NCO_outpQ;
//     float rds_normBandwidth = 0.003;
//     float rds_phaseAdjust = 0.0;
//     float rds_ncoScale = 0.5;


//     std::vector<float> dem_mixer;
//     std::vector<float> dem_rds_LPF_coeffs;
//     std::vector<float> dem_rds_LPF_filtered;
//     std::vector<float> dem_state_rds_LPF(values.num_Taps-1, 0.0);


//     float dem_resamplerFs = 2375 * values.SPS;
//     float dem_resamplerFc = std::min((values.audio_expan / values.audio_decim) *(2375.0/2), (2375.0/2));
// 	std::vector<float> dem_rds_resamp_coeffs(values.num_Taps * values.audio_expan, 0.0);
//     std::vector<float> dem_rds_resamp_filtered;
//     std::vector<float> dem_state_rds_resamp(values.num_Taps-1, 0.0);

// 	std::vector<float> RRC_coeffs;
// 	std::vector<float> RF_h;
// 	// std::vector<float> final_coeffs;

// 	std::vector<float> i_data, q_data;
// 	std::vector<float> filt_i, filt_q;
// 	std::vector<float> demod;
// 	std::vector<float> state_i(values.num_Taps-1, 0.0);
// 	std::vector<float> state_q(values.num_Taps-1, 0.0);
// 	// std::vector<float> state_mono(num_Taps-1, 0.0);

// 	// std::cerr<<"TEST"<<std::endl;
// 	float prev_i = 0.0;
// 	float prev_q = 0.0;
    
//     while (!done) {
//         if (demod_ptr) {
//             // process rds data
//         }
//         // rds synchronization, extraction, etc.....
//         // after processing, notify the condition_variable to possibly (just maybe) wake up the rf thread ????
//         BPFCoeffs(RDSFb, RDSFe, values.IF_Fs, values.num_Taps, rds_BPF_coeffs); // last is output
//         //BPF CARRIER RECOVERY
//         BPFCoeffs(CR_RDSFb, CR_RDSFe, values.IF_Fs, values.num_Taps, CR_rds_BPF_coeffs); // last is output
//         impulseResponseLPF(values.IF_Fs, 3000, values.num_Taps, dem_rds_LPF_coeffs); // last is output
//         gainimpulseResponseLPF(dem_resamplerFs, dem_resamplerFc, values.num_Taps, dem_rds_resamp_coeffs, values.rds_up); // 4th is output

//         // DEMOD DATA
//         std::shared_ptr<std::vector<float>> demod = tsQueue.wait_and_pop(); // Block until data is available
        
//         conv_ds_fast(rds_filtered, demod, rds_BPF_coeffs, 1, state_rds_bb);
//         // std::cerr<<"after rds 0"<<std::endl;
//         rds_nonlin.resize(rds_filtered.size());
//         for(int i = 0; i < rds_filtered.size(); i++) {
//             rds_nonlin[i] = rds_filtered[i] * rds_filtered[i];
//         }
//         // std::cerr<<"after non linearity"<<std::endl;
//         //ALL PASS FILTER
//         delayBlock(rds_filtered, rds_processed_delay, rds_state_delay);
//         conv_ds_fast(CR_rds_filtered, rds_nonlin, CR_rds_BPF_coeffs, 1, CR_state_rds);
//         // std::cerr<<"tst"<<std::endl;
//         //PLL
//         fmPll(CR_rds_filtered, rds_NCO_outp, rds_NCO_outpQ, rds_lockInFreq, values.IF_Fs, rds_ncoScale, rds_phaseAdjust, rds_normBandwidth, rds_state); 		
//         // std::cerr<<"after pll"<<std::endl;

//         dem_mixer.resize(CR_rds_filtered.size());
//         for(int i = 0; i<CR_rds_filtered.size();i++){
//             dem_mixer[i] = 3*rds_NCO_outp[i] * rds_processed_delay[i];
//         }
        
//         conv_ds_fast(dem_rds_LPF_filtered, dem_mixer, dem_rds_LPF_coeffs, 1, dem_state_rds_LPF);
//         conv_rs(dem_rds_resamp_filtered, dem_rds_LPF_filtered, dem_rds_resamp_coeffs, values.rds_down, values.rds_up, dem_state_rds_LPF);

//         impulseResponseRootRaisedCosine(RRC_coeffs, rds_fs_rat_res, values.num_Taps);
//         conv_ds_fast(outp_rrc, dem_rds_resamp_filtered, RRC_coeffs, 1, state_rrc);
//         std::vector<float> bits;
//         // std::cerr<<"before get bits"<<std::endl;
//         get_bits(outp_rrc, values.SPS, bits);
//         // std::cerr<<"after getting bits"<<std::endl;
    
//     }
// }


