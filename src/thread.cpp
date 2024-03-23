#include "dy4.h"
#include "thread.h"
#include "mode.h"
#include "filter.h"
#include "RFfront.h"

std::atomic<bool> done{false};

template <typename T>
class ThreadSafeQueue {
private:
    mutable std::mutex m;                   // Mutex to protect access to the queue
    std::queue<std::shared_ptr<T>> q;       // Standard queue wrapped inside the thread-safe queue
    std::condition_variable cv;             // Condition variable for notifying waiting threads

public:
    // Enqueues an element by copying it into a shared_ptr and adding it to the queue
    void push(T value) {
        std::shared_ptr<T> data(std::make_shared<T>(std::move(value))); // Create a shared_ptr to the new data
        std::lock_guard<std::mutex> lock(m);  // Lock the mutex during the queue operation
        q.push(data);                         // Push the data onto the queue
        cv.notify_one();                      // Notify one waiting thread that there is new data available
    }

    // Waits for an element to be available and pops it from the queue
    std::shared_ptr<T> wait_and_pop() {
        std::unique_lock<std::mutex> lock(m); // Unique lock allows the lock to be temporarily released
        cv.wait(lock, [this]{ return !q.empty(); }); // Wait until the queue is not empty
        std::shared_ptr<T> result = q.front(); // Get the front element
        q.pop();                              // Remove the element from the queue
        return result;                         // Return the data to the caller
    }

    // Tries to pop an element from the queue without waiting
    bool try_pop(T& value) {
        std::lock_guard<std::mutex> lock(m);  // Lock the mutex during the queue operation
        if(q.empty()) {
            return false;                     // Return false if the queue is empty
        }
        value = std::move(*q.front());        // Move the value from the front of the queue into the provided variable
        q.pop();                              // Remove the element from the queue
        return true;                          // Return true to indicate a value was popped
    }

    // Checks if the queue is empty
    bool empty() const {
        std::lock_guard<std::mutex> lock(m);  // Lock the mutex to synchronize access to the empty check
        return q.empty();                     // Return whether the queue is empty
    }

    size_t size() const {
        std::lock_guard<std::mutex> lock(m);
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

ThreadSafeQueue<std::vector<float>> tsQueue; // Global instance of the thread-safe queue

// // Function representing the work of the RF thread (the producer)
void rf_thread(int mode)  {                        // Continue producing until done is true
    Mode values;
    values.configMode(mode);
    std::vector<float> i_data, q_data;
	std::vector<float> filt_i, filt_q;
    std::vector<float> state_i(values.num_Taps, 0.0);
	std::vector<float> state_q(values.num_Taps, 0.0);
    std::vector<float> RF_h;
    std::vector<float> demod;
    float prev_i = 0.0;
    float prev_q = 0.0;
    gainimpulseResponseLPF(values.RF_Fs, values.RF_Fc, values.num_Taps, RF_h, values.audio_expan);
    while (true){
        for(unsigned int block_id = 0; ; block_id++){
            //std::cerr<<"Block id "<<block_id<<std::endl;
            std::vector<float> block_data(values.BLOCK_SIZE);
            readStdinBlockData(values.BLOCK_SIZE, block_id, block_data);
            if((std::cin.rdstate()) != 0){
                std::cerr<<"End of input stream reached" << std::endl;
                tsQueue.print_contents();
                exit(1);
            }
            std::cerr<<"Block id "<<block_id<<std::endl;
            split_audio_iq(block_data, i_data, q_data);
            conv_ds_fast(filt_i, i_data, RF_h, values.rf_decim, state_i);
            conv_ds_fast(filt_q, q_data, RF_h, values.rf_decim, state_q);
            FM_demod(filt_i, filt_q, prev_i, prev_q, demod);
            // std::cerr<<"size before pushing: "<<tsQueue.size()<<std::endl;
            tsQueue.push(demod);
        }
    }
        // Push the produced data onto the queue

}

// // Function representing the work of the audio thread (the consumer)
// void audio_thread() {
//     while (!done) {                           // Continue consuming until done is true
//         std::shared_ptr<std::vector<float>> data_ptr = tsQueue.wait_and_pop(); // Wait for and pop data from the queue
//         consume_data(*data_ptr);              // Consume data (placeholder function)
//     }
// }

// // Entry point of the program
// int main() {
//     std::thread producer(rf_thread);          // Create the producer thread
//     std::thread consumer(audio_thread);       // Create the consumer thread

//     // Other code to control the threads, for example to set 'done' when needed

//     producer.join();                          // Wait for the producer thread to finish
//     consumer.join();                          // Wait for the consumer thread to finish

//     return 0;
// }

// // Placeholder function for data production
// bool done = false; // Global flag to control the thread loop execution

std::vector<float> produce_data(int mode) {
    Mode values;
    values.configMode(mode);
    std::vector<float> i_data, q_data;
	std::vector<float> filt_i, filt_q;
    std::vector<float> state_i(values.num_Taps, 0.0);
	std::vector<float> state_q(values.num_Taps, 0.0);
    std::vector<float> RF_h;
    std::vector<float> demod;
    float prev_i = 0.0;
    float prev_q = 0.0;
    gainimpulseResponseLPF(values.RF_Fs, values.RF_Fc, values.num_Taps, RF_h, values.audio_expan);
    while (true){
        for(unsigned int block_id = 0; ; block_id++){
            std::cerr<<"Block id "<<block_id<<std::endl;
            std::vector<float> block_data(values.BLOCK_SIZE);
            readStdinBlockData(values.BLOCK_SIZE, block_id, block_data);
            if((std::cin.rdstate()) != 0){
                std::cerr<<"End of input stream reached" << std::endl;
                exit(1);
            }
            split_audio_iq(block_data, i_data, q_data);
            conv_ds_fast(filt_i, i_data, RF_h, values.rf_decim, state_i);
            conv_ds_fast(filt_q, q_data, RF_h, values.rf_decim, state_q);
            FM_demod(filt_i, filt_q, prev_i, prev_q, demod);
            
        }
    }
}

// Placeholder function for data consumption
void consume_data(const std::vector<float>& data) {
    // Your data consumption logic
}
