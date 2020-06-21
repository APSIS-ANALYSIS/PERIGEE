"""
Ingrid Lan. June 2020.

# Generates the 'inflow_fourier_series.txt' file required by PERIGEE. Fourier smoothing is performed in the same way
# as in SimVascular, where the time series in the *.flow file is first linearly interpolated to NUM_BCT_PTS. The
# interpolated time series is then transformed into Fourier domain and subsequently transformed back into time with
# NUM_FOURIER_MODES.

"""

import glob
import matplotlib.pyplot as plt
import numpy as np
from scipy.fftpack import fft


# PERIGEE simulation directory
pg_sim_dir = '/home/ingridxlan/Documents/Ingrid/Solvers/PERIGEE_testing/PigQ/1.2M_pg/'


# Parameters in svpre for inflow interpolation & Fourier smoothing
NUM_BCT_PTS = 201
NUM_FOURIER_MODES = 10

# ================================================================================================================== #

"""
Linear interplation & FFT to generate coefficients Q_coeff, a, b in the following two forms
Q_recon1 = sum[ Q_coeff * exp(i n w t) ]
Q_recon2 = sum[a * cos(n w t) + b * sin(n w t)]

"""
def fourier_coeff(t_data, q_data, period, N_interp):

    # ==== SV starts the interpolation from the first time point in the inflow.flow file ====
    # This shifts the time waveform to start from 0
    t_interp = np.linspace(t_data[0], period, N_interp)
    q_interp = np.interp(t_interp, t_data, q_data, period=period)

    N = len(t_interp)
    if N % 2 != 0:      # odd
        num_fourier_modes = 1 + (N - 1) / 2
        print(N)
    else:
        num_fourier_modes = 1 + N / 2

    fft_inflow = fft(q_interp)
    Q_coeff = fft_inflow[ : num_fourier_modes] / N
    Q_coeff[1 : ] *= 2.0

    a = np.real(Q_coeff)
    b = -np.imag(Q_coeff)

    return (t_interp, Q_coeff, a, b)


"""
Fourier reconstruction in 2 different forms that should yield the same waveform.
Q_recon1 = sum[ Q_coeff * exp(i n w t) ]
Q_recon2 = sum[a * cos(n w t) + b * sin(n w t)]

"""
def fourier_recon(t_interp, Q_coeff, a, b, period, N_fourier_modes):

    Q_recon1 = np.zeros(len(t_interp), dtype=np.complex128)
    Q_recon2 = np.zeros(len(t_interp), dtype=np.complex128)

    for n in range(N_fourier_modes):
        w = 2.0 * np.pi / period
        Q_recon1 += Q_coeff[n] * np.exp(complex(0.0, 1.0) * n * w * t_interp)
        Q_recon2 += a[n] * np.cos(n * w * t_interp) + b[n] * np.sin(n * w * t_interp)

    return (w, Q_recon1, Q_recon2)

# ================================================================================================ #


if __name__ == "__main__":

    flow_file = glob.glob(pg_sim_dir + '*.flow')[0]
    q_read = np.loadtxt(flow_file)
    t_data =  q_read[:, 0]
    q_data = -q_read[:, 1]
    period =  t_data[-1]

    t_interp, Q_coeff, a, b = fourier_coeff(t_data, q_data, period, NUM_BCT_PTS)
    w, Q_recon1, Q_recon2   = fourier_recon(t_interp, Q_coeff, a, b, period, NUM_FOURIER_MODES)

    # Plot the raw data & two reconstructions (identical)
    plt.figure()
    plt.plot(t_data, q_data, 'ro', label='Raw Data')
    plt.plot(t_interp, Q_recon1, 'k--', label='Recon 1')
    plt.plot(t_interp, Q_recon2, 'b:', label='Recon 2')

    plt.ylabel('Q (mL/s)')
    plt.xlabel('Time (s)')
    plt.savefig(pg_sim_dir + 'inflow_fourier_smoothing.png')
    plt.close()

    # Write out 'inflow_fourier_series.txt' for perigee
    with open(pg_sim_dir + 'inflow_fourier_series.txt', 'w') as outfile:
        outfile.write('# num_fourier_modes - 1 / fundamental frequency (w) / period\n')
        outfile.write(str(NUM_FOURIER_MODES - 1) + ' ' + str(w) + ' ' + str(period) + '\n')
        
        outfile.write('\n# Coefficients a_n with length = num_fourier_modes\n')
        a_n_str = ''
        for n in range(NUM_FOURIER_MODES):
            a_n_str += str(a[n]) + ' '
        a_n_str += '\n'
        outfile.write(a_n_str)

        outfile.write('\n# Coefficients b_n with length = num_fourier_modes\n')
        b_n_str = ''
        for n in range(NUM_FOURIER_MODES):
            b_n_str += str(b[n]) + ' '
        b_n_str += '\n'
        outfile.write(b_n_str)