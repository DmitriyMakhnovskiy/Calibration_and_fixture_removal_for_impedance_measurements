#
# Function for unwrapping the phase sweep of S21(f) or S12(f) obtained from the 3-term error model.
# The sign of the unwrapped phase slope as well as the delay time along the cable are calculated.
#

import numpy as np
import math
from scipy.optimize import minimize
from scipy import stats

def unwrap(f, Re, Im):
    # f - frequency array, Hz
    # Re - Re[S21(f)] array
    # Im - Im[S21(f)] array
    angle = list(map(math.atan2, Im, Re))  # original phase array with jumps
    unwangle = list(map(math.atan2, Im, Re))  # initial values for the unwrapped phase array
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
    
    # Additional vertical displacements will be required to compensate the horizontal shift of the phase pieces
    # after the jumps. Since these displacements depend on slopes of the phase pieces,
    # we will call them "gradient shifts".
    grad = []  # array of the gradient shifts
    slopex = []  # frequency points between the left and right jumps
    slopey = []  # phase values between the left and right jumps
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
    num = 0
    denum = 0
    for i in range(anglenum):
        num = num + f[i] * unwangle[i]
        denum = denum + f[i]**2
    ratio = num / denum
    slopesign = np.sign(ratio)  # sign of the slope of the unwrapped phase
    dt = np.abs(ratio) / (2.0 * np.pi)  # delay time calculated from the slope of the unwrapped phase

    # This function is used for the fine adjustment of the delay time by minimization of the frequency dispersion
    # profile of |Im(S21(f)/exp(-i*w*dt))| by dt. The minimisation functional is the sum of |Im(S21(f)/exp(-i*w*dt))|.
    # over all frequency points.
    def optfun(x):
        amp = []  # amplitude array of the imaginary part
        for i in range(anglenum):
            value = np.abs(Re[i] * np.sin(2.0 * np.pi * f[i] * x) + Im[i] * np.cos(2.0 * np.pi * f[i] * x))
            amp.append(value)
        value = 0.0
        for i in range(anglenum):
            value = value + amp[i]
        return value

    # Minimization of optfun by Nelder-Mead simplex method.
    # https://en.wikipedia.org/wiki/Nelder-Mead_method
    res = minimize(optfun, dt, method='nelder-mead', options={'xatol': 1e-20})
    dt = res.x[0]  # final value of the delay time obtained after the optimisation

    return slopesign, dt, unwangle