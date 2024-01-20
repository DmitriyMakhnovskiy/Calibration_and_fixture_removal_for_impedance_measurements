#
# Cubic spline of the initial data (array CalS) and recalculating it on a new scale (Stimulus)
#

import numpy as np
from scipy.interpolate import CubicSpline

def cubspl(CalS, Stimulus):
    N = len(Stimulus)
    freq = CalS[:, 0]  # Frequencies at which the self-made terminations were measured
    Re = CalS[:, 1]  # Real part
    Im = CalS[:, 2]  # Imaginary part
    ReCub = CubicSpline(freq, Re)  # Cubic spline interpolation
    ImCub = CubicSpline(freq, Im)  # Cubic spline interpolation
    ReNew = ReCub(Stimulus)  # New values over the actual frequency sweep (Stimulus)
    ImNew = ImCub(Stimulus)  # New values over the actual frequency sweep (Stimulus)
    ReNew.reshape([N, ], order='F')  # Column-major (Fortran-style) order in memory
    ImNew.reshape([N, ], order='F')  # Column-major (Fortran-style) order in memory
    S = np.empty([N, ], dtype=np.complex128, order='F')  # Column of complex numbers
    S.real = ReNew
    S.imag = ImNew
    return S