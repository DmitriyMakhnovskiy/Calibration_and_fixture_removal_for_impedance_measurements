#
# Delay time algorithm
#
# The algorithm calculates the delay time along a fixture/sample under test (wire or microstip above the ground plane,
# coaxial probe, waveguide structure, or free space) using the unwrapped phase of the S11 or S21 parameter
# measured on VNA.
#
# Updated 24.01.2024
# DYK team: http://dykteam.com/
#

import numpy as np
from Acquisition import acq
from PhaseUnwrapping import unwrap
import matplotlib.pyplot as plt

print('')
print('*******************************************< Program operation >**********************************************')
print('')
print('(1) VNA together with the cables must be already calibrated with a coaxial calibration kit.')
print('     For a wide frequency range, use SOLT calibration (TOSM in terms of Rohde&Schwarz).')
print('')
print('(2) Activate the VNA calibration (and the 2-port deembedding if required).')
print('')
print('(3) Measure the dispersion of S11 or S21 for the fixture/sample under test and save it as S.CSV.')
print('    The file must have three columns: frequency (Hz), Re[S], Im[S]')
print('')
print('(4) First, the program will plot the initial phase dispersion where you have to choose the frequency range ')
print('    with the most linear behavior. It may be with or without jumps. If with jumps, ensure that all phase lines')
print('    are approximately parallel. If it is not possible, choose the frequency range without jumps.')
print('')
print('(5) The delay time will be calculated from the unwrapped phase withing the selected frequency range.')
print('    Output files: Phase_initial.csv and Phase_unwrapped.csv')
print('')
print('**************************************************************************************************************')
print('')

folder = input("Enter/paste the path of your folder where the file S.CSV is located: ")

deliminput = input('Please choose the delimiter used in S.CSV ( , ; tab space): ')

print('')

Stimulus0, MS0 = acq('S.CSV', folder, deliminput)

MS0 = np.array(MS0)
Stimulus0 = np.array(Stimulus0)
N0 = len(Stimulus0)

Re0 = MS0.real
Im0 = MS0.imag
Ph0 = np.angle(MS0)

phaseoutput = np.column_stack((Stimulus0, Ph0))  # 2D array
np.savetxt(folder + '\\' + 'Phase_initial.CSV', phaseoutput, delimiter=',')
print('The initial phase has been saved in Phase_initial.csv')
print('')

input("""Now the program will draw the initial phase dispersion before unwrapping.
Note at which value, pi/2 (~1.5) or pi (~3), the phase jumps.
Later, choose the phase factor 1 if the phase jumps at pi/2, and the phase factor 2 if the phase jumps at pi.
Push Enter when finish reading.""")
print('')

font1 = {'family': 'serif', 'color': 'black', 'weight': 'bold', 'size': 15}
font2 = {'family': 'serif', 'color': 'black', 'weight': 'normal', 'size': 10}
plt.plot(Stimulus0, Ph0, 'k')
plt.title('Phase dispersion', fontdict=font1)
plt.xlabel('Frequency, Hz', fontdict=font2)
plt.ylabel('Phase, radians', fontdict=font2)
plt.show()

phase_factor = int(input('Please enter the phase factor (1 or 2), as explained above: '))

flag = 0
while flag == 0:
    fstart = np.double(input('Enter the start frequency in Hz: '))
    fstop = np.double(input('Enter the stop frequency in Hz: '))

    if Stimulus0[0] <= fstart <= fstop <= Stimulus0[N0 - 1]:
        flag = 1
    else:
        print('')
        print('Oppps... You entered a wrong sequence of frequencies')
        print('It must be: fmin <= fstart <= fstop <= fmax')
        print('Please start again')
        print('')

# New frequency range (narrower) that includes only the linear behavior of phase
Stimulus = []
Re = []
Im = []
MS = []
for i in range(0, N0 - 1):
    if fstart <= Stimulus0[i] <= fstop:
        Stimulus.append(Stimulus0[i])
        Re.append(Re0[i])
        Im.append(Im0[i])
        MS.append(MS0[i])

N = len(Stimulus)
Stimulus = np.array(Stimulus)
Re = np.array(Re)
Im = np.array(Im)
MS = np.array(MS)
Ph = np.angle(MS)

font1 = {'family': 'serif', 'color': 'black', 'weight': 'bold', 'size': 15}
font2 = {'family': 'serif', 'color': 'black', 'weight': 'normal', 'size': 10}
plt.plot(Stimulus, Ph, 'k')
plt.title('Phase dispersion', fontdict=font1)
plt.xlabel('Frequency, Hz', fontdict=font2)
plt.ylabel('Phase, radians', fontdict=font2)
plt.show()

slopesign, dt, unwrapped_phase, initial_phase = unwrap(Stimulus, Re, Im, phase_factor)  # phase unwrapping
unwrappedphaseoutput = np.column_stack((Stimulus, unwrapped_phase))  # 2D array
np.savetxt(folder + '\\' + 'Phase_unwrapped.csv', unwrappedphaseoutput, delimiter=',')

font1 = {'family': 'serif', 'color': 'black', 'weight': 'bold', 'size': 15}
font2 = {'family': 'serif', 'color': 'black', 'weight': 'normal', 'size': 10}
plt.plot(Stimulus, unwrapped_phase, 'k')
plt.title('Unwrapped phase', fontdict=font1)
plt.xlabel('Frequency, Hz', fontdict=font2)
plt.ylabel('Phase, radians', fontdict=font2)
plt.show()

print('The unwrapped phase has been saved in Phase_unwrapped.csv')
print('Delay time = ', dt / 1.0e-12, 'ps')
print('')