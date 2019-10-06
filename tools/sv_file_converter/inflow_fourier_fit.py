# $ sudo pip install symfit


from symfit import parameters, variables, sin, cos, Fit
import numpy as np
import matplotlib.pyplot as plt

def fourier_series(x, f, n=0):
    """
    Returns a symbolic fourier series of order `n`.

    :param n: Order of the fourier series.
    :param x: Independent variable
    :param f: Frequency of the fourier series
    """
    # Make the parameter objects for all the terms
    cos_a = parameters(','.join(['a{}'.format(i) for i in range(0, n + 1)]))
    a0 = cos_a[0]
    cos_a = cos_a[1 : ]
    sin_b = parameters(','.join(['b{}'.format(i) for i in range(1, n + 1)]))
    # Construct the series
    series = a0 + sum(ai * cos(i * f * x) + bi * sin(i * f * x)
                     for i, (ai, bi) in enumerate(zip(cos_a, sin_b), start=1))
    return series

# ================================================================================================ #

num_fourier_modes = 6
x, y = variables('x, y')
w, = parameters('w')
model_dict = {y: fourier_series(x, f=w, n=num_fourier_modes)}
print(model_dict)


q_read = np.loadtxt('inflow_scaleFactor_0.8.flow')
t_data =  q_read[:, 0]
q_data = -q_read[:, 1]
period =  t_data[-1]


fit = Fit(model_dict, x=t_data, y=q_data)
fit_result = fit.execute()
print(fit_result)

model_params = fit_result.model.params

# Plot the result
plt.plot(t_data, q_data, 'ro')
plt.plot(t_data, fit.model(x=t_data, **fit_result.params).y, 'k--')
plt.ylabel('Q (mL/s)')
plt.xlabel('Time (s)')
plt.savefig('inflow_fourier_fit.png')

# Write out 'inflow_fourier_series.txt'
with open('inflow_fourier_series.txt', 'w') as outfile:
    outfile.write('# num_fourier_modes / fundamental frequency (w) / period\n')
    w = fit_result.value(model_params[-1]) # last parameter
    outfile.write(str(num_fourier_modes) + ' ' + str(w) + ' ' + str(period) + '\n')
    
    outfile.write('\n# Coefficients a_n with length = num_fourier_modes + 1\n')
    a_n_str = ''
    for i in range(num_fourier_modes + 1):
        a_n_str += str(fit_result.value(model_params[i])) + ' '
    a_n_str += '\n'
    outfile.write(a_n_str)

    outfile.write('\n# Coefficients b_n with length = num_fourier_modes + 1\n')
    b_n_str = '0.0 '
    for i in range(num_fourier_modes + 1, len(model_params) - 1):
        b_n_str += str(fit_result.value(model_params[i])) + ' '
    b_n_str += '\n'
    outfile.write(b_n_str)

