#
# Acquisition of the real and imaginary parts from a file and creating an array of complex numbers
#

import numpy as np

def acq(filename, folder, delim):
    address = folder + '\\' + filename  #  full address of the file
    data = np.genfromtxt(address, delimiter=delim)  # reading the csv file
    Freq = data[:, 0]  # first column - frequency
    N = len(Freq)  # number of frequency points
    Real = data[:, 1]  # second column - real part
    Imag = data[:, 2]  # third column - imaginary part
    Freq.reshape([N, ], order='F')  # Column-major (Fortran-style) order in memory
    Real.reshape([N, ], order='F')  # Column-major (Fortran-style) order in memory
    Imag.reshape([N, ], order='F')  # Column-major (Fortran-style) order in memory
    S = np.empty([N, ], dtype=np.complex128, order='F')  # Column of complex numbers
    S.real = Real
    S.imag = Imag
    return Freq, S
