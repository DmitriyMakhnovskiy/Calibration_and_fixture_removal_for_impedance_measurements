#
# Probe deembedding algorithm
#
# The algorithm creates complete S-parameter models of the probes (one or two) based on measurements
# of the reflection coefficients from the self-made SHORT, OPEN and LOAD terminations attached to the probes free end.
# If required (old or simplified VNA), the probe(s) will be deembedded for the 1-port or 2-port measurements.
#
# For the 1-port measurements, the device under test is connected to the port 1. Only S11 is measured.
# With the created S2P file of the probe, it can be deembedded automatically on a modern VNA, or numerically.
#
# For the 2-port measurements, the device under test is connected in series between the two probes.
# All S-parameters are measured.
# With the created S2P files of the probes, they can be deembedded automatically on a modern VNA, or numerically.
#
# 20 January 2024
# DYK team: http://dykteam.com/

print('')
print('****************************************< Program operation >**********************************************')
print('')
print('(1) VNA together with the cables must be already calibrated with a coaxial calibration kit.')
print('     For a wide frequency range, use SOLT calibration (TOSM in terms of Rohde&Schwarz)')
print('')
print('(2) Activate the VNA calibration. You should see a sign indicating that the calibration is enabled.')
print('')
print('(3) If using non-ideal self-made SOL terminations, you must provide their S11/22 dispersions:')
print('     S11S.csv - dispersion of the reflection from the SHORT termination connected to the probe 1')
print('     S11O.csv - dispersion of the reflection from the OPEN termination connected to the probe 1')
print('     S11L.csv - dispersion of the reflection from the LOAD termination connected to the probe 1')
print('     S22S.csv - dispersion of the reflection from the SHORT termination connected to the probe 2')
print('     S22O.csv - dispersion of the reflection from the OPEN termination connected to the probe 2')
print('     S22L.csv - dispersion of the reflection from the LOAD termination connected to the probe 2')
print('     Each file has three columns: frequency (Hz), real part, imaginary part')
print('')
print('(4) Measure the following reflection coefficients for the 3-term model and save them in the following files:')
print('     MS11S.csv - reflection measured from the probe 1 terminated with SHORT')
print('     MS11O.csv - reflection measured from the probe 1 terminated with OPEN')
print('     MS11L.csv - reflection measured from the probe 1 terminated with LOAD')
print('     MS22S.csv - reflection measured from the probe 2 terminated with SHORT')
print('     MS22O.csv - reflection measured from the probe 2 terminated with OPEN')
print('     MS22L.csv - reflection measured from the probe 2 terminated with LOAD')
print('     Each file has three columns: frequency (Hz), real part, imaginary part')
print('')
print(' (5) (Only for old VNA) Measure S-parameters from the device under test and save them in the following files:')
print('      MS11.csv - reflection measured from the first port (1-port device or the first port of a 2-port device)')
print('      MS22.csv - reflection measured from the second port of a 2-port device under test')
print('      MS21.csv - transmission measured through a 2-port device under test')
print('      MS12.csv - transmission measured through a 2-port device under test')
print('      Each file has three columns: frequency (Hz), real part, imaginary part')
print('      Note: These measurements are not required if using a modern VNA with S2P deembedding files')
print('')
print(' (6) Main program outputs:')
print('     (a) Reflection (S11A) or transmission (S21A) from/through the device under test after deembedding')
print('         will be calculated and saved in the files S11A.csv and S21A.csv')
print('         Files format: freq (Hz), Re[S], Im[S], |S|, log10(|S|)')
print('     (b) ZA_from_S11A.csv or ZA_from_S21A.csv - actual impedance of the device under test calculated')
print('         Files format: freq (Hz), Re[ZA], Im[ZA], |ZA|, log10(freq), log10(|ZA|)')
print('     (c) Probe_1.s2p and Probe_2.s2p - complete S-parameter models of the probes saved in the files')
print('         These files can be uploaded on VNA for automatic deembedding')
print('************************************************************************************************************')
print('')
exec(open('Deembedding.py').read())  # Deembedding the probe