#

# Comp Eng 3DY4 (Computer Systems Integration Project)

#

# Copyright by Nicola Nicolici

# Department of Electrical and Computer Engineering

# McMaster University

# Ontario, Canada

#



import numpy as np

import math

def fmPll(pllIn, freq, Fs, ncoScale=1.0, phaseAdjust=0.0, normBandwidth=0.01, state=None):
    # Previous state initialization
    if state is None:
        state = {
            'integrator': 0.0,
            'phaseEst': 0.0,
            'feedbackI': 1.0,
            'feedbackQ': 0.0,
            'trigOffset': 0
        }

    Cp = 2.666
    Ci = 3.555

    Kp = (normBandwidth) * Cp
    Ki = (normBandwidth * normBandwidth) * Ci

    ncoOut = np.empty(len(pllIn) + 1)
    ncoOut[0] = 1.0

    integrator = state['integrator']
    phaseEst = state['phaseEst']
    feedbackI = state['feedbackI']
    feedbackQ = state['feedbackQ']
    trigOffset = state['trigOffset']
    prev_ncoOut = state['prev_ncoOut']


    ncoOut[0] = prev_ncoOut
    
    for k in range(len(pllIn)):
        errorI = pllIn[k] * (+feedbackI)
        errorQ = pllIn[k] * (-feedbackQ)
        errorD = math.atan2(errorQ, errorI)

        integrator = integrator + Ki * errorD
        phaseEst = phaseEst + Kp * errorD + integrator

        trigOffset += 1
        trigArg = 2 * math.pi * (freq / Fs) * trigOffset + phaseEst
        feedbackI = math.cos(trigArg)
        feedbackQ = math.sin(trigArg)
        ncoOut[k + 1] = math.cos(trigArg * ncoScale + phaseAdjust)
        if(k < 7 or k > len(pllIn) - 7):
            print("index: ", k, "\ttrigArg: ", round(trigArg, 6), "\tfeedbackI: ", round(feedbackI,6),"\tncoOut at index", k,  round(ncoOut[k],6),  "\tfeedbackQ: ", round(feedbackQ,6), "\ttrigOffset", round(trigOffset,6), "\tPLLin at k: ", round(pllIn[k],6),"\tphaseEst: ", round(phaseEst,6), "integrator", round(integrator, 6))
    prev_ncoOut = ncoOut[-1]
    state = {
        'integrator': integrator,
        'phaseEst': phaseEst,
        'feedbackI': feedbackI,
        'feedbackQ': feedbackQ,
        'trigOffset': trigOffset,
        'prev_ncoOut':prev_ncoOut
    }

    return ncoOut, state

# def fmPll(pllIn, freq, Fs, ncoScale = 1.0, phaseAdjust = 0.0, normBandwidth = 0.01, state=None):
    """
    pllIn          array of floats
                   input signal to the PLL (assume known frequency)

    freq           float
                   reference frequency to which the PLL locks

    Fs             float
                   sampling rate for the input/output signals

    ncoScale       float
                   frequency scale factor for the NCO output

    phaseAdjust    float
                   phase adjust to be added to the NCO output only

    normBandwidth  float
                   normalized bandwidth for the loop filter
                   (relative to the sampling rate)

    state          dictionary
                   dictionary to store and update internal state

    """
    # Scale factors for proportional/integrator terms
    # These scale factors were derived assuming a damping factor of 0.707 (1 over square root of 2)
    # There is no oscillator gain and no phase detector gain
    Cp = 2.666
    Ci = 3.555

    # Gain for the proportional term
    Kp = normBandwidth * Cp
    # Gain for the integrator term
    Ki = normBandwidth * normBandwidth * Ci

    # Initialize internal state if not provided
    if state is None:
        state = {
            'integrator': 0.0,
            'phaseEst': 0.0,
            'feedbackI': 1.0,
            'feedbackQ': 0.0,
            'ncoOut': np.empty(len(pllIn) + 1),
            'trigOffset': 0
        }
        state['ncoOut'][0] = 1.0

    # Process each sample
    for k in range(len(pllIn)):
        # Phase detector
        errorI = pllIn[k] * (+state['feedbackI'])  # Complex conjugate of the input
        errorQ = pllIn[k] * (-state['feedbackQ'])  # Feedback complex exponential

        # Four-quadrant arctangent discriminator for phase error detection
        errorD = math.atan2(errorQ, errorI)

        # Loop filter
        state['integrator'] = state['integrator'] + Ki * errorD

        # Update phase estimate
        state['phaseEst'] = state['phaseEst'] + Kp * errorD + state['integrator']

        # Internal oscillator
        state['trigOffset'] += 1
        trigArg = 2 * math.pi * (freq / Fs) * state['trigOffset'] + state['phaseEst']
        state['feedbackI'] = math.cos(trigArg)
        state['feedbackQ'] = math.sin(trigArg)
        state['ncoOut'][k + 1] = math.cos(trigArg * ncoScale + phaseAdjust)
        if(k < 7 or k > len(pllIn) - 7):
            print("index: ", k, "\ttrigArg: ", round(trigArg, 6), "\tfeedbackI: ", round(state['feedbackI'],6),"\tncoOut at index", k,  round(state['ncoOut'][k],6),  "\tfeedbackQ: ", round(state['feedbackQ'],6), "\ttrigOffset", round(state['trigOffset'],6), "\tPLLin at k: ", round(pllIn[k],6),"\tphaseEst: ", round(state['phaseEst'],6), "integrator", round(state['integrator'], 6))

            #print("feedbackQ: ", state['feedbackQ'], "\tintegrator", state['integrator'], "\tphaseEst", state['phaseEst'])
    # Update the last element of ncoOut
    state['ncoOut'][-1] = math.cos((2 * math.pi * (freq / Fs) * (state['trigOffset'] + 1)) * ncoScale + phaseAdjust)
    print('\n')
    # Return output and updated state
    return state['ncoOut'], state
# def fmPll(pllIn, freq, Fs, ncoScale = 1.0, phaseAdjust = 0.0, normBandwidth = 0.01):
# 	"""

# 	pllIn 	 		array of floats

# 					input signal to the PLL (assume known frequency)

# 	freq 			float

# 					reference frequency to which the PLL locks

# 	Fs  			float

# 					sampling rate for the input/output signals

# 	ncoScale		float

# 					frequency scale factor for the NCO output

# 	phaseAdjust		float

# 					phase adjust to be added to the NCO output only

# 	normBandwidth	float

# 					normalized bandwidth for the loop filter

# 					(relative to the sampling rate)

# 	state 			to be added


# 	"""

# 	# scale factors for proportional/integrator terms

# 	# these scale factors were derived assuming the following:

# 	# damping factor of 0.707 (1 over square root of 2)

# 	# there is no oscillator gain and no phase detector gain

# 	Cp = 2.666
# 	Ci = 3.555

# 	# gain for the proportional term
# 	Kp = (normBandwidth)*Cp
# 	# gain for the integrator term
# 	Ki = (normBandwidth*normBandwidth)*Ci

# 	# output array for the NCO
# 	ncoOut = np.empty(len(pllIn)+1)

# 	# initialize internal state
# 	integrator = 0.0
# 	phaseEst = 0.0
# 	feedbackI = 1.0
# 	feedbackQ = 0.0
# 	ncoOut[0] = 1.0
# 	trigOffset = 0

# 	# note: state saving will be needed for block processing



# 	for k in range(len(pllIn)):

# 		# phase detector
# 		errorI = pllIn[k] * (+feedbackI)  # complex conjugate of the

# 		errorQ = pllIn[k] * (-feedbackQ)  # feedback complex exponential



# 		# four-quadrant arctangent discriminator for phase error detection
# 		errorD = math.atan2(errorQ, errorI)

# 		# loop filter
# 		integrator = integrator + Ki*errorD

# 		# update phase estimate
# 		phaseEst = phaseEst + Kp*errorD + integrator

# 		# internal oscillator
# 		trigOffset += 1
# 		trigArg = 2*math.pi*(freq/Fs)*(trigOffset) + phaseEst
# 		feedbackI = math.cos(trigArg)
# 		feedbackQ = math.sin(trigArg)
# 		ncoOut[k+1] = math.cos(trigArg*ncoScale + phaseAdjust)

# 	# for stereo only the in-phase NCO component should be returned
# 	# for block processing you should also return the state
# 	return ncoOut,
	# for RDS add also the quadrature NCO component to the output



if __name__ == "__main__":



	pass
