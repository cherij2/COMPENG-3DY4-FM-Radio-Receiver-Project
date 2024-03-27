#
# Comp Eng 3DY4 (Computer Systems Integration Project)
#
# Copyright by Nicola Nicolici(da goat)
# Department of Electrical and Computer Engineering
# McMaster University
# Ontario, Canada
#

import numpy as np
import math

def manchester_encoding(rrc_after, bits, temp, temp2):

    for i in range(rrc_after-1):
        #hard coding state saving cuz im dumb
        if i == 0:
            if (temp > rrc_after[0] and temp > temp2):
                bits.append(1)  # local max
                temp2 = temp
            if (temp < rrc_after[0] and temp < temp2):
                bits.append(0)  # local min
                temp2 = temp
        if i == 1:
            if (rrc_after[0]  > rrc_after[1] and rrc_after[0]  > temp2):
                bits.append(1)  # local max
                temp2 = temp
            if (rrc_after[0]  < rrc_after[1] and rrc_after[0]  < temp2):
                bits.append(0)  # local min
                temp = rrc_after[i]
        if i == len(rrc_after) - 1:
            temp = rrc_after[i] 
            temp2 = rrc_after[i-1] 
        if (temp > rrc_after[i + 1] and temp > rrc_after[i - 1]):
            bits.append(1)  # local max
            temp = rrc_after[i]
        if (temp < rrc_after[i + 1] and temp < rrc_after[i - 1]):
            bits.append(0)  # local min
            temp = rrc_after[i]


    #state saving

    return bits
#fs = 71250 samples/sec, 30 samples/symbol
#1 symbol is half of a period of sinusoid
# check every 30 samples if we still positive
def manchester_encoding_fast(rrc_after, bits, SPS, temp):
    temp = 0
    
    for i in range(rrc_after-1):
        if(rrc_after[i] < temp):
            






def manchester_decode(bits):
    manch_decoded = []
    state = -3  
    for i in range(0, len(bits) - 1, 2):  
        if bits[i] == 0 and bits[i + 1] == 1:
            manch_decoded.append(0)
        elif bits[i] == 1 and bits[i + 1] == 0:
            manch_decoded.append(1)
        elif (bits[i] == 1 and bits[i + 1] == 1) or (bits[i] == 0 and bits[i + 1] == 0):
            continue
        
        
        # State saving
        if i == len(bits) - 1:
            state = bits[i]  
        elif i == 0 and state != -3:
            if state == 0 and bits[0] == 1:
                manch_decoded.append(0)
            elif state == 1 and bits[0] == 0:
                manch_decoded.append(1)
    return manch_decoded
