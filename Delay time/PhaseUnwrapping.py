#
# Function for unwrapping the phase sweep of S21(f) or S11(f)
#

import numpy as np
from scipy import stats

def unwrap(f, Re, Im, phase_factor):
    # f - frequency array, Hz
    # Re - Re[S21(f)] array
    # Im - Im[S21(f)] array
    angle = np.arctan2(Im, Re)  # original phase array with jumps
    unwangle = angle  # initial values for the unwrapped phase array
    anglenum = len(angle)  # number of the phase values that must coincide with the number of frequencies

    # A phase jump consists only of two points - "left" and "right", which have opposite signs. No other points between.
    left_jump = []  # array of the indexes for the jump "start points"
    right_jump = []  # array of the indexes for the jump "stop points"

    # Searching for the phase jumps.
    for i in range(anglenum - 2):
        sign = np.sign(angle[i] * angle[i + 1])
        if sign < 0 and np.abs(angle[i + 2]) <= np.abs(angle[i + 1]):  # criteria for selecting a jump
            left_jump.append(i)
            right_jump.append(i + 1)

    jumpnum = len(right_jump)  # number of the jump points
    print('')
    print('Number of phase jumps = ', jumpnum)
    print('')

    for i in range(jumpnum):
        for j in range(anglenum):
            if right_jump[i] <= j:
                unwangle[j] = unwangle[j] - phase_factor * np.pi

     #  Calculating the slope of the unwrapped phase
    ratio = stats.linregress(f, unwangle)[0]
    slopesign = np.sign(ratio)  # sign of the slope of the unwrapped phase
    dt = np.abs(ratio) / (2.0 * np.pi)  # delay time (s) calculated from the slope of the unwrapped phase

    return slopesign, dt, unwangle, angle