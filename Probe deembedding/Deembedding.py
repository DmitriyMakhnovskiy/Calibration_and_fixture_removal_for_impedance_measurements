#
# Deembedding one or two probes
# If two probes, their termination SOL files must be measured on the same set of frequency points
#

import numpy as np
import math
from Acquisition import acq
from scipy import linalg
from PhaseUnwrapping import unwrap
from CubSpline import cubspl
from Jumps import jumps

# Initial configuration
flag = 0
while flag == 0:
    probe_flag = int(input('Will you use one or two probes? Enter 1 or 2 respectively: '))
    print('')
    if probe_flag != 1 and probe_flag != 2:
        flag = 0
    else:
        flag = 1

folder = input("Enter/paste the path of your folder where all files will be located: ")
print('')

deliminput = input('Please choose the delimiter used in the input files ( , ; tab space): ')
print('')

delimoutput = input('''
VNA (e.g. Rohde&Schwarz) may prefer the tab delimiter in S2P file. Please check when uploading the file on VNA. 
Please choose the delimiter used in the output files (, ; tab space ): ''')
print('')

flag = 0
while flag == 0:
    ideal = int(input('Will you use ideal termination model? If yes, enter 1, otherwise 0: '))
    print('')
    if ideal != 1 and ideal != 0:
        flag = 0
    else:
        flag = 1

flag = 0
while flag == 0:
    VNA = int(input('Will you use a modern VNA with automatic S2P file deembedding? If yes, enter 1, otherwise 0: '))
    print('')
    if VNA != 1 and VNA != 0:
        flag = 0
    else:
        flag = 1

print('')

#  Measurements when the probe is terminated with SHORT (S), OPEN (O), and LOAD (L)
Stimulus, MS11S = acq('MS11S.CSV', folder, deliminput)
_, MS11O = acq('MS11O.CSV', folder, deliminput)
_, MS11L = acq('MS11L.CSV', folder, deliminput)

if probe_flag == 2:
    _, MS22S = acq('MS22S.CSV', folder, deliminput)
    _, MS22O = acq('MS22O.CSV', folder, deliminput)
    _, MS22L = acq('MS22L.CSV', folder, deliminput)

if ideal == 0:  # non-ideal terminations

# The reflections for the SOL self-made standards could be measured at frequencies other than the actual sweep.
# That is why these reflections from the SOL standards must be recalculated over the actual frequency sweep points.
# To do so, we use cubic spline interpolation. Stimulus - actual frequency points
    # Non-ideal terminations for the probe 1
    address = folder + '\\' + 'S11S.CSV'
    CalS11S = np.genfromtxt(address, delimiter=deliminput)
    S11S = cubspl(CalS11S, Stimulus)  # Complex values recalculated over the actual frequency sweep points.

    address = folder + '\\' + 'S11O.CSV'
    CalS11O = np.genfromtxt(address, delimiter=deliminput)
    S11O = cubspl(CalS11O, Stimulus)  # Complex values recalculated over the actual frequency sweep points.

    address = folder + '\\' + 'S11L.CSV'
    CalS11L = np.genfromtxt(address, delimiter=deliminput)
    S11L = cubspl(CalS11L, Stimulus)  # Complex values recalculated over the actual frequency sweep points.

    if probe_flag == 2:  # non-ideal terminations for the probe 2
        address = folder + '\\' + 'S22S.CSV'
        CalS22S = np.genfromtxt(address, delimiter=deliminput)
        S22S = cubspl(CalS22S, Stimulus)  # Complex values recalculated over the actual frequency sweep points.

        address = folder + '\\' + 'S22O.CSV'
        CalS22O = np.genfromtxt(address, delimiter=deliminput)
        S22O = cubspl(CalS22O, Stimulus)  # Complex values recalculated over the actual frequency sweep points.

        address = folder + '\\' + 'S22L.CSV'
        CalS22L = np.genfromtxt(address, delimiter=deliminput)
        S22L = cubspl(CalS22L, Stimulus)  # Complex values recalculated over the actual frequency sweep points.

else:  # ideal terminations (option ideal = 1)
    S11S = []
    S11O = []
    S11L = []
    for i in range(len(Stimulus)):
        S11S.append(complex(-1.0, 0.0))
        S11O.append(complex(1.0, 0.0))
        S11L.append(complex(0.0, 0.0))

    if probe_flag == 2:
        S22S = []
        S22O = []
        S22L = []
        for i in range(len(Stimulus)):
            S22S.append(complex(-1.0, 0.0))
            S22O.append(complex(1.0, 0.0))
            S22L.append(complex(0.0, 0.0))

# Function for solving the 3-term error model
# AS - pre-saved reflection coefficient of the SHORT termination
# AO - pre-saved reflection coefficient of the OPEN termination
# AL - pre-saved reflection coefficient of the LOAD termination
# AMS - measured reflection coefficient of the probe when it is terminated with SHORT
# AMO - measured reflection coefficient of the probe when it is terminated with OPEN
# AML - measured reflection coefficient of the probe when it is terminated with LOAD
def three_term_error(AS, AO, AL, AMS, AMO, AML):
    A = [[1.0, AS * AMS, AS],  # matrix of the 3-term model
         [1.0, AO * AMO, AO],
         [1.0, AL * AML, AL]]
    b = [AMS, AMO, AML]
    vector = linalg.solve(A, b)  # Parameters of the network: S11, S22, z = S21**2 - S11 * S22 (S21 = S12)
    return vector

# Arrays S11(f) = x and S22(f) = y
# S21 = S12 = +/-(z + x * y)**0.5; the sign has to be chosen from the condition of the negative phase slope
x = []
y = []
z = []
for i in range(len(Stimulus)):  # Loop over the actual frequency points
    AS = S11S[i]  # SHORT dispersion
    AO = S11O[i]  # OPEN dispersion
    AL = S11L[i]  # LOAD dispersion
    AMS = MS11S[i]  # Probe terminated with SHORT
    AMO = MS11O[i]  # Probe terminated with OPEN
    AML = MS11L[i]  # Probe terminated with LOAD
    xyz = three_term_error(AS, AO, AL, AMS, AMO, AML)  # Solution (vector) of the 3-term error model
    x.append(xyz[0])  # Dispersion of x
    y.append(xyz[1])  # Dispersion of y
    z.append(xyz[2])  # Dispersion of z

x = np.array(x)
y = np.array(y)
z = np.array(z)

# S21 and S12 will be calculated as transvar^0.5 after unwrapping the phase of the complex transvar
transvar = z + x * y  # definition of transvar (transmission variable)
realtransvar = transvar.real  # real part
imagtransvar = transvar.imag  # imaginary part
phase = np.arctan2(imagtransvar, realtransvar)  # atan2 function at https://en.wikipedia.org/wiki/Atan2

phaseoutput = np.column_stack((Stimulus, phase))  # 2D array
np.savetxt(folder + '\\' + 'Probe_1_phase_initial.CSV', phaseoutput, delimiter=delimoutput)
print('The probe 1 initial phase has been saved in Probe_1_phase_initial.csv')
print('')

input("""Now open Probe_1_phase_initial.CSV and check at which value, pi/2 (~1.5) or pi (~3), the phase jumps.
Later, choose the phase factor 1 if the phase jumps at pi/2, and the phase factor 2 if the phase jumps at pi.
Push Enter when finish reading.""")
print('')

phase_factor = int(input("Please enter the phase factor (1 or 2) as explained above: "))
print('')

number_jumps = jumps(realtransvar, imagtransvar)
if number_jumps == 0:
    slopesign = np.sign(phase[int(len(Stimulus) / 2.0)])
else:
    slopesign, dt, unwrapped_phase, initial_phase = unwrap(Stimulus, realtransvar, imagtransvar, phase_factor)
    unwrappedphaseoutput = np.column_stack((Stimulus, unwrapped_phase))  # 2D array
    np.savetxt(folder + '\\' + 'Probe_1_phase_unwrapped.CSV', unwrappedphaseoutput, delimiter=delimoutput)
    print('The probe 1 unwrapped phase has been saved in Probe_1_phase_unwrapped.csv')
    print('Delay time along the probe 1 = ', dt / 1.0e-12, 'ps')
    print('')

# Choosing the proper sign of the slope: it must be negative
if slopesign >= 0:
    unwrapped_phase = -unwrapped_phase

# Calculating S21 and S12 from the complex transvar = z + x * y with the unwrapped phase
# Trans - transmission
realTrans = []
imagTrans = []
for i in range(len(Stimulus)):
    realTrans.append(((np.absolute(transvar[i]))**0.5) * np.cos(unwrapped_phase[i] / 2.0))
    imagTrans.append(((np.absolute(transvar[i]))**0.5) * np.sin(unwrapped_phase[i] / 2.0))

realTrans = np.array(realTrans)
imagTrans = np.array(imagTrans)

# S2P matrix of the probe 1
#   Header: # Hz S RI R 50.0'
#   9 columns: frequency, Re[S11], Im[S11], Re[S21], Im[S21], Re[S12], Im[S12], Re[S22], Im[S22].
# For a probe, S21 = S12.
Probe1_S2P = np.column_stack((Stimulus, x.real, x.imag, realTrans, imagTrans, realTrans, imagTrans, y.real, y.imag))

# Saving S2P file on PC
np.savetxt(folder + '\\' + 'Probe_1.S2P', Probe1_S2P, delimiter=delimoutput, header='Hz S RI R 50.00')
print('The probe 1 S2P file has been saved in Probe_1.s2p')
print('')

if probe_flag == 1 and VNA == 0:
    #  Forming outputs in the case of a simplified or old VNA (option 0) for a one port measurement
    S11A = []
    ZA = []
    LinMagS11A = []
    LogMagS11A = []
    LinMagZA = []
    LogMagZA = []
    LogMagStimulus = []
    Stimulus, MS11 = acq('MS11.CSV', folder, deliminput)
    Stimulus = np.array(Stimulus)  # array format
    MS11 = np.array(MS11)  # array format
    for i in range(len(Stimulus)):

        LogMagStimulusvalue = math.log10(abs(Stimulus[i]))

        S11Avalue = MS11[i] - x[i]
        S11Avalue = S11Avalue / (z[i] + x[i] * y[i] + y[i] * (MS11[i] - x[i]))
        LinMagS11Avalue = abs(S11Avalue)
        LogMagS11Avalue = math.log10(LinMagS11Avalue)

        ZAvalue = 50.0 * (complex(1.0, 0.0) + S11Avalue) / (complex(1.0, 0.0) - S11Avalue)
        LinMagZAvalue = abs(ZAvalue)
        LogMagZAvalue = math.log10(LinMagZAvalue)

        S11A.append(S11Avalue)
        ZA.append(ZAvalue)
        LinMagS11A.append(LinMagS11Avalue)
        LogMagS11A.append(LogMagS11Avalue)
        LinMagZA.append(LinMagZAvalue)
        LogMagZA.append(LogMagZAvalue)
        LogMagStimulus.append(LogMagStimulusvalue)

    S11A = np.array(S11A)  # array format for calculating the real and imaginary parts
    ZA = np.array(ZA)  # array format for calculating the real and imaginary parts

    OutputS11A = np.column_stack((Stimulus, S11A.real, S11A.imag, LinMagS11A, LogMagS11A))
    np.savetxt(folder + '\\' + 'S11A.CSV', OutputS11A, delimiter=delimoutput)
    print('The actual S11A of the device under test has been saved in S11A.csv')
    print('')

    OutputZA = np.column_stack((Stimulus, ZA.real, ZA.imag, LinMagZA, LogMagStimulus, LogMagZA))
    np.savetxt(folder + '\\' + 'ZA_from_S11A.CSV', OutputZA, delimiter=delimoutput)
    print('The actual ZA of the device under test calculated from S11A has been saved in ZA_from_S11A.csv')
    print('')

if probe_flag == 2:
    x = []
    y = []
    z = []
    for i in range(len(Stimulus)):  # Loop over the actual frequency points
        AS = S22S[i]  # SHORT dispersion
        AO = S22O[i]  # OPEN dispersion
        AL = S22L[i]  # LOAD dispersion
        AMS = MS22S[i]  # Probe terminated with SHORT
        AMO = MS22O[i]  # Probe terminated with OPEN
        AML = MS22L[i]  # Probe terminated with LOAD
        xyz = three_term_error(AS, AO, AL, AMS, AMO, AML)  # Solution (vector) of the 3-term error model
        x.append(xyz[0])  # Dispersion of x
        y.append(xyz[1])  # Dispersion of y
        z.append(xyz[2])  # Dispersion of z

    x = np.array(x)
    y = np.array(y)
    z = np.array(z)

    # S21 and S12 will be calculated as transvar^0.5 after unwrapping the phase of the complex transvar
    transvar = z + x * y  # definition of transvar
    realtransvar = transvar.real
    imagtransvar = transvar.imag
    phase = np.arctan2(imagtransvar, realtransvar)  # atan2 function at https://en.wikipedia.org/wiki/Atan2

    phaseoutput = np.column_stack((Stimulus, phase))  # 2D array
    np.savetxt(folder + '\\' + 'Probe_2_phase_initial.CSV', phaseoutput, delimiter=delimoutput)
    print('The probe 2 initial phase has been saved in Probe_2_phase_initial.csv')
    print('')

    input("""Now open Probe_2_phase_initial.CSV and check at which value, pi/2 (~1.5) or pi (~3), the phase jumps.
    Later, choose the phase factor 1 if the phase jumps at pi/2, and the phase factor 2 if the phase jumps at pi.
    Push Enter when finish reading.""")
    print('')

    phase_factor = int(input("Please enter the phase factor (1 or 2) as explained above: "))
    print('')

    number_jumps = jumps(realtransvar, imagtransvar)
    if number_jumps == 0:
        slopesign = np.sign(phase[int(len(Stimulus) / 2.0)])
    else:
        slopesign, dt, unwrapped_phase, initial_phase = unwrap(Stimulus, realtransvar, imagtransvar, phase_factor)
        unwrappedphaseoutput = np.column_stack((Stimulus, unwrapped_phase))  # 2D array
        np.savetxt(folder + '\\' + 'Probe_2_phase_unwrapped.CSV', unwrappedphaseoutput, delimiter=delimoutput)
        print('The probe 2 unwrapped phase has been saved in Probe_2_phase_unwrapped.csv')
        print('Delay time along the probe 2 = ', dt / 1.0e-12, 'ps')
        print('')

    # Choosing the proper sign of the slope: it must be negative
    if slopesign >= 0:
        unwrapped_phase = -unwrapped_phase

    # Calculating S21 and S12 from the complex transvar = z + x * y with the unwrapped phase
    # Trans - transmission
    realTrans = []
    imagTrans = []
    for i in range(len(Stimulus)):
        realTrans.append(((np.absolute(transvar[i])) ** 0.5) * np.cos(unwrapped_phase[i] / 2.0))
        imagTrans.append(((np.absolute(transvar[i])) ** 0.5) * np.sin(unwrapped_phase[i] / 2.0))

    realTrans = np.array(realTrans)
    imagTrans = np.array(imagTrans)

    # S2P matrix of the probe 2
    #   Header: # Hz S RI R 50.0'
    #   9 columns: frequency, Re[S11], Im[S11], Re[S21], Im[S21], Re[S12], Im[S12], Re[S22], Im[S22].
    # For a probe, S21 = S12.
    Probe2_S2P = np.column_stack((Stimulus, x.real, x.imag, realTrans, imagTrans, realTrans, imagTrans, y.real, y.imag))

    # Saving S2P file on PC
    np.savetxt(folder + '\\' + 'Probe_2.S2P', Probe2_S2P, delimiter=delimoutput, header='Hz S RI R 50.00')
    print('The probe 2 S2P file has been saved in Probe_2.s2p')
    print('')

if probe_flag == 2 and VNA == 0:
    # S-parameters of the device under test
    Stimulus, MS11 = acq('MS11.CSV', folder, deliminput)
    _, MS22 = acq('MS22.CSV', folder, deliminput)
    _, MS21 = acq('MS21.CSV', folder, deliminput)
    _, MS12 = acq('MS12.CSV', folder, deliminput)
    MS11 = np.array(MS11)
    MS22 = np.array(MS22)
    MS21 = np.array(MS21)
    MS12 = np.array(MS12)
    Stimulus = np.array(Stimulus)

    # S-parameter model of the first probe
    a = Probe1_S2P[:,1]
    b = Probe1_S2P[:,2]
    a.reshape([len(Stimulus), ], order='F')  # Column-major (Fortran-style) order in memory
    b.reshape([len(Stimulus), ], order='F')  # Column-major (Fortran-style) order in memory
    S11_1 = np.empty([len(Stimulus), ], dtype=np.complex128, order='F')  # Column of complex numbers
    S11_1.real = a
    S11_1.imag = b
    S11_1 = np.array(S11_1)

    a = Probe1_S2P[:,3]
    b = Probe1_S2P[:,4]
    a.reshape([len(Stimulus), ], order='F')  # Column-major (Fortran-style) order in memory
    b.reshape([len(Stimulus), ], order='F')  # Column-major (Fortran-style) order in memory
    S21_1 = np.empty([len(Stimulus), ], dtype=np.complex128, order='F')  # Column of complex numbers
    S21_1.real = a
    S21_1.imag = b
    S21_1 = np.array(S21_1)

    a = Probe1_S2P[:,5]
    b = Probe1_S2P[:,6]
    a.reshape([len(Stimulus), ], order='F')  # Column-major (Fortran-style) order in memory
    b.reshape([len(Stimulus), ], order='F')  # Column-major (Fortran-style) order in memory
    S12_1 = np.empty([len(Stimulus), ], dtype=np.complex128, order='F')  # Column of complex numbers
    S12_1.real = a
    S12_1.imag = b
    S12_1 = np.array(S12_1)

    a = Probe1_S2P[:,7]
    b = Probe1_S2P[:,8]
    a.reshape([len(Stimulus), ], order='F')  # Column-major (Fortran-style) order in memory
    b.reshape([len(Stimulus), ], order='F')  # Column-major (Fortran-style) order in memory
    S22_1 = np.empty([len(Stimulus), ], dtype=np.complex128, order='F')  # Column of complex numbers
    S22_1.real = a
    S22_1.imag = b
    S22_1 = np.array(S22_1)

    # S-parameter model of the second probe
    a = Probe2_S2P[:,1]
    b = Probe2_S2P[:,2]
    a.reshape([len(Stimulus), ], order='F')  # Column-major (Fortran-style) order in memory
    b.reshape([len(Stimulus), ], order='F')  # Column-major (Fortran-style) order in memory
    S11_2 = np.empty([len(Stimulus), ], dtype=np.complex128, order='F')  # Column of complex numbers
    S11_2.real = a
    S11_2.imag = b
    S11_2 = np.array(S11_2)

    a = Probe2_S2P[:,3]
    b = Probe2_S2P[:,4]
    a.reshape([len(Stimulus), ], order='F')  # Column-major (Fortran-style) order in memory
    b.reshape([len(Stimulus), ], order='F')  # Column-major (Fortran-style) order in memory
    S21_2 = np.empty([len(Stimulus), ], dtype=np.complex128, order='F')  # Column of complex numbers
    S21_2.real = a
    S21_2.imag = b
    S21_2 = np.array(S21_2)

    a = Probe2_S2P[:,5]
    b = Probe2_S2P[:,6]
    a.reshape([len(Stimulus), ], order='F')  # Column-major (Fortran-style) order in memory
    b.reshape([len(Stimulus), ], order='F')  # Column-major (Fortran-style) order in memory
    S12_2 = np.empty([len(Stimulus), ], dtype=np.complex128, order='F')  # Column of complex numbers
    S12_2.real = a
    S12_2.imag = b
    S12_2 = np.array(S12_2)

    a = Probe2_S2P[:,7]
    b = Probe2_S2P[:,8]
    a.reshape([len(Stimulus), ], order='F')  # Column-major (Fortran-style) order in memory
    b.reshape([len(Stimulus), ], order='F')  # Column-major (Fortran-style) order in memory
    S22_2 = np.empty([len(Stimulus), ], dtype=np.complex128, order='F')  # Column of complex numbers
    S22_2.real = a
    S22_2.imag = b
    S22_2 = np.array(S22_2)

    S21A = []
    ZA = []
    LinMagS21A = []
    LogMagS21A = []
    LinMagZA = []
    LogMagZA = []
    LogMagStimulus = []
    for i in range(len(Stimulus)):
        PM = np.array([[MS12[i] * MS21[i] - MS11[i] * MS22[i], MS22[i]],
                        [-MS11[i], complex(1.0, 0.0)]] / MS12[i])
        InvP1 = np.array([[complex(1.0, 0.0), -S22_1[i]],
                        [S11_1[i], S21_1[i] * S12_1[i] - S11_1[i] * S22_1[i]]] / S21_1[i])
        InvP2 = np.array([[complex(1.0, 0.0), -S11_2[i]],
                        [S22_2[i], S12_2[i] * S21_2[i] - S11_2[i] * S22_2[i]]] / S12_2[i])
        PA = InvP2 @ (PM @ InvP1)

        S21Avalue = PA[0][0] - PA[0][1] * PA[1][0] / PA[1][1]
        ZAvalue = 100.0 * (complex(1.0, 0.0) - S21Avalue) / S21Avalue
        S21A.append(S21Avalue)
        ZA.append(ZAvalue)

        LinMagS21Avalue = abs(S21Avalue)
        LogMagS21Avalue = math.log10(LinMagS21Avalue)
        LinMagZAvalue = abs(ZAvalue)
        LogMagZAvalue = math.log10(LinMagZAvalue)
        LogMagStimulusvalue = math.log10(abs(Stimulus[i]))

        LinMagS21A.append(LinMagS21Avalue)
        LogMagS21A.append(LogMagS21Avalue)
        LinMagZA.append(LinMagZAvalue)
        LogMagZA.append(LogMagZAvalue)
        LogMagStimulus.append(LogMagStimulusvalue)

    S21A = np.array(S21A)  # array format for calculating the real and imaginary parts
    ZA = np.array(ZA)  # array format for calculating the real and imaginary parts

    OutputS21A = np.column_stack((Stimulus, S21A.real, S21A.imag, LinMagS21A, LogMagS21A))
    np.savetxt(folder + '\\' + 'S21A.CSV', OutputS21A, delimiter=delimoutput)
    print('The actual transmission S21A of DUT has been saved in S21A.csv')
    print('')

    OutputZA = np.column_stack((Stimulus, ZA.real, ZA.imag, LinMagZA, LogMagStimulus, LogMagZA))
    np.savetxt(folder + '\\' + 'ZA_from_S21A.CSV', OutputZA, delimiter=delimoutput)
    print('The actual impedance ZA of DUT calculated from S21A has been saved in ZA_from_S21A.csv')
    print('')