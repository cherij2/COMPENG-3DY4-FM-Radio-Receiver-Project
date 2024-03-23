#ifndef DY4_THREAD_H
#define DY4_THREAD_H

#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <math.h>
#include <mutex>
#include <condition_variable>
#include <memory>
#include <queue>
#include <thread>
#include <iostream>
#include <vector>
#include <atomic>

template<typename T>
class ThreadSafeQueue;
std::vector<float> produce_data();
void consume_data();


#endif