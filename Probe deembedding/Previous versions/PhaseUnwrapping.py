#
# Function for unwrapping the phase sweep of S21(f) or S12(f) obtained from the 3-term error model.
# The sign of the unwrapped phase slope as well as the delay time along the cable are calculated.
#

import numpy as np
import math
from scipy import stats

def unwrap(f, Re, Im):
    # f - frequency array, Hz
    # Re - Re[S21(f)] array
    # Im - Im[S21(f)] array
    angle = list(map(math.atan2, Im, Re))  # original phase array with jumps
    unwangle = angle  # initial values for the unwrapped phase array
    anglenum = len(angle)  # number of the phase values that must coincides with the number of frequencies

    # A phase jump consists only of two points - "left" and "right", which have opposite signs. No other points between.
    left_jump = []  # array of the indexes of the jump "start points"
    right_jump = []  # array of the indexes of the jump "stop points"
    yshift = []  # array of the vertical shifts used in the unwrapping method

    # Searching for the phase jumps.
    for i in range(anglenum - 2):
        sign = np.sign(angle[i] * angle[i + 1])
        if sign < 0 and np.abs(angle[i + 2]) <= np.abs(angle[i + 1]):  # criteria for the selecting a jump
            left_jump.append(i)
            right_jump.append(i + 1)

    jumpnum = len(right_jump)  # number of the jump points
    print('')
    print('Number of phase jumps = ', jumpnum)
    print('')

    # Additional vertical displacements will be required to compensate the horizontal shift of the phase pieces
    # after the jumps. Since these displacements depend on slopes of the phase pieces,
    # we will call them "gradient shifts".
    grad = []  # array of the gradient shifts
    slopex = []  # frequency points between the left and right jumps
    slopey = []  # phase values between the left and right jumps
    if jumpnum != 0:
        if left_jump[0] != 0:  # phase sweep starts not with a jump
            for i in range(left_jump[0] + 1):
                slopex.append(f[i])
                slopey.append(angle[i])
            a = stats.linregress(slopex, slopey)[0]  # least-square for calculating the slope
            b = stats.linregress(slopex, slopey)[1]  # least-square for calculating the interception
            grad.append(a * f[right_jump[0]] + b - angle[left_jump[0]])  # first element has been created
        else:  # if the phase sweep start with a jump
            for i in range(right_jump[0], left_jump[1] + 1):
                slopex.append(f[i])
                slopey.append(angle[i])
            a = stats.linregress(slopex, slopey)[0]
            b = stats.linregress(slopex, slopey)[1]
            grad.append(a * f[right_jump[0]] + b - angle[left_jump[0]])

    if jumpnum != 0:
        for i in range(1, jumpnum):  # the rest of gradient shifts
            slopex = []
            slopey = []
            for j in range(right_jump[i - 1], left_jump[i]):
                slopex.append(f[j])
                slopey.append(angle[j])
            a = stats.linregress(slopex, slopey)[0]
            b = stats.linregress(slopex, slopey)[1]
            grad.append(a * f[right_jump[i]] + b - angle[left_jump[i]])

        for i in range(jumpnum):
            yshift.append(angle[right_jump[i]] - angle[left_jump[i]] - grad[i])  # total vertical shifts

        # Phase unwrapping: shifting the phase values down or up to unwrap the values into a straight line
        for i in range(jumpnum):
            for j in range(right_jump[i], anglenum):
                unwangle[j] = unwangle[j] - yshift[i]

    #  Calculating the slope of the unwrapped phase
    ratio = stats.linregress(f, unwangle)[0]
    slopesign = np.sign(ratio)  # sign of the slope of the unwrapped phase
    dt = np.abs(ratio) / (2.0 * np.pi)  # delay time calculated from the slope of the unwrapped phase

    return slopesign, dt, unwangle