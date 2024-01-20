#
# Function counts the number of phase jumps
#

import numpy as np

def jumps(Re, Im):
    # f - frequency array, Hz
    # Re - Re[S21(f)] array
    # Im - Im[S21(f)] array
    angle = np.arctan2(Im, Re)  # original phase array with jumps
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

    return len(right_jump)  # number of the jump points