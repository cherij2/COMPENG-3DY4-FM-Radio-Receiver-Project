// #include "dy4.h"
// #include "thread.h"

// std::atomic<bool> done{false};

// template <typename T>
// class ThreadSafeQueue {
// private:
//     mutable std::mutex m;                   // Mutex to protect access to the queue
//     std::queue<std::shared_ptr<T>> q;       // Standard queue wrapped inside the thread-safe queue
//     std::condition_variable cv;             // Condition variable for notifying waiting threads

// public:
//     // Enqueues an element by copying it into a shared_ptr and adding it to the queue
//     void push(T value) {
//         std::shared_ptr<T> data(std::make_shared<T>(std::move(value))); // Create a shared_ptr to the new data
//         std::lock_guard<std::mutex> lock(m);  // Lock the mutex during the queue operation
//         q.push(data);                         // Push the data onto the queue
//         cv.notify_one();                      // Notify one waiting thread that there is new data available
//     }

//     // Waits for an element to be available and pops it from the queue
//     std::shared_ptr<T> wait_and_pop() {
//         std::unique_lock<std::mutex> lock(m); // Unique lock allows the lock to be temporarily released
//         cv.wait(lock, [this]{ return !q.empty(); }); // Wait until the queue is not empty
//         std::shared_ptr<T> result = q.front(); // Get the front element
//         q.pop();                              // Remove the element from the queue
//         return result;                         // Return the data to the caller
//     }

//     // Tries to pop an element from the queue without waiting
//     bool try_pop(T& value) {
//         std::lock_guard<std::mutex> lock(m);  // Lock the mutex during the queue operation
//         if(q.empty()) {
//             return false;                     // Return false if the queue is empty
//         }
//         value = std::move(*q.front());        // Move the value from the front of the queue into the provided variable
//         q.pop();                              // Remove the element from the queue
//         return true;                          // Return true to indicate a value was popped
//     }

//     // Checks if the queue is empty
//     bool empty() const {
//         std::lock_guard<std::mutex> lock(m);  // Lock the mutex to synchronize access to the empty check
//         return q.empty();                     // Return whether the queue is empty
//     }
// };

// // ==================================

// ThreadSafeQueue<std::vector<float>> tsQueue; // Global instance of the thread-safe queue

// // Function representing the work of the RF thread (the producer)
// void rf_thread() {
//     while (!done) {                           // Continue producing until done is true
//         std::vector<float> fm_demodulated_data = produce_data(); // Produce data (placeholder function)
//         tsQueue.push(std::move(fm_demodulated_data)); // Push the produced data onto the queue
//     }
// }

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

// std::vector<float> produce_data(int mode) {
//     std::vector<float> i_data, q_data;
// 	std::vector<float> filt_i, filt_q;
//     std::vector<float> demod;
//     std::vector<float> state_i(num_Taps, 0.0);
// 	std::vector<float> state_q(num_Taps, 0.0);
//     while (true){
//         for(unsigned int block_id = 0; ; block_id++){
//             std::vector<float> block_data(BLOCK_SIZE)
//             readStdinBlockData(BLOCK_SIZE, block_id, block_data);
//             if((std::cin.rdstate()) != 0){
//                 std::cerr<<"End of input stream reached" << std::endl;
//                 exit(1);
//             }
//             split_audio_iq(block_data, i_data, q_data);
//             conv_ds_fast(filt_i, i_data, RF_h, rf_decim, state_i);
//             conv_ds_fast(filt_q, q_data, RF_h, rf_decim, state_q);
//             FM_demod(filt_i, filt_q, prev_i, prev_q, demod);
//         return demod;
//         }
//     }
    

// }

// // Placeholder function for data consumption
// void consume_data(const std::vector<float>& data) {
//     // Your data consumption logic
// }
