#
# Calculation of the impedance dispersion from the measured and corrected S11 or S21
# The delay time along the sample will be required to correct S11/S22.
#
# Yujie Zhao, University of St. Andrews, Scotland, 18 March 2022
# Dmitriy Makhnovskiy, Sensing Materials Technology Ltd, UK
# DYK team: http://dykteam.com/
#

import numpy as np
import matplotlib.pyplot as plt

print('')
print('Your S-parameter (S11 or S21) must be saved as S.csv with three columns: freq (Hz), Re, Im')
print('Remove any headers in the file.')
print('')

folder = input("Enter/paste the path of your folder where S.csv is located: ")
deliminput = input('Please choose the delimiter used in the input files ( , ; tab space): ')
parameter = int(input("Will you calculate the impedance from S11 or S21? Enter 11 or 21 respectively: "))

address = folder + '\\' + 'S.csv'  #  full address of the file
data = np.genfromtxt(address, delimiter=deliminput)  # reading the csv file
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

Freq = np.array(Freq)

while True:
    delta = float(input('Enter the delay time along the sample in ps: '))
    print('')
    delta = delta * 1.0e-12
    SS = []
    for i in range(N):
        SS.append(S[i] * np.exp(complex(0.0, 1.0) * 2.0 * np.pi * Freq[i] * delta))  # corrected S-parameter

    SS = np.array(SS)
    if parameter == 11:
        Z = 50.0 * (complex(1.0, 0.0) + SS) / (complex(1.0, 0.0) - SS)
    elif parameter == 21:
        Z = 100.0 * (complex(1.0, 0.0) - SS) / SS
    else:
        print('Did you use S11 or S21? Please run the program again.')
        exit()

    Z = np.array(Z)

    S_corrected = np.column_stack((Freq, SS.real, SS.imag, np.angle(SS)))
    Z_corrected = np.column_stack((Freq, Z.real, Z.imag, np.abs(Z)))
    np.savetxt(folder + '\\' + 'S_corrected.CSV', S_corrected, delimiter=',')
    np.savetxt(folder + '\\' + 'Z_corrected.CSV', Z_corrected, delimiter=',')

    font1 = {'family': 'serif', 'color': 'black', 'weight': 'bold', 'size': 15}
    font2 = {'family': 'serif', 'color': 'black', 'weight': 'normal', 'size': 10}
    plt.plot(Freq, Z.real, 'g')
    plt.plot(Freq, Z.imag, 'r')
    plt.title('Impedance dispersion: Re - green, Im - red', fontdict=font1)
    plt.xlabel('Frequency, Hz', fontdict=font2)
    plt.ylabel('Z, Ohms', fontdict=font2)
    plt.grid()
    plt.show()