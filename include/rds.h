#ifndef DY4_RDS_H
#define DY4_RDS_H

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

void impulseResponseRootRaisedCosine(std::vector<float> &h, float Fs, int N_taps);


void get_bits (const std::vector<float> &rrc_after, int sps, std::vector<float> &bits);

void preDataProcessing(int mode);


#endif