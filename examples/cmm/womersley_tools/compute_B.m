function [B_n, Q_n, G_n] = compute_B(flow, rho, nu_s, c_n, gamma_n, g_n, T, R, n_modes)

% base angular frequency
omega = 2 * pi / T;
assignin('base', 'omega', omega);

% Time vector
N = floor(length(flow));
t = linspace(0, T, N);

% Discrete fast Fourier transform
Q_n = fft(flow);

% Compute the Fourier coefficients for non-negative frequencies
Q_n = Q_n / N;
Q_n = Q_n(1 : N / 2 + 1);
Q_n(2 : end - 1) = 2 * Q_n(2 : end - 1);

% Truncate to first n modes (including the steady 0th mode)
Q_n = Q_n(1 : n_modes);

% ------- Plot the input vs. reconstructed flow --------
figure; hold on; plot(t, flow, 'b-o');

flow_recon = zeros(1, length(t));
for k = 1 : n_modes
    n = k - 1;
    flow_recon = flow_recon + Q_n(k) * exp(1j * n * omega * t);
end

plot(t, real(flow_recon), 'r-');
xlabel('Time (s)'); ylabel('Flow (mL/s)')
legend('Input', 'Reconstructed');
set(gca, 'FontSize', 12)

% Elasticity factor. Only entries starting at index 2 are meaningful
G_n = ( 2 + gamma_n * (2 * nu_s - 1) ) ./ ( gamma_n .* (2 * nu_s - g_n) );

% Solve for constant B needed to define the oscillatory flow
% Only entries starting at index 2 are meaningful
B_n = Q_n * rho .* c_n ./ ( pi * R^2 * (1 - G_n .* g_n) );
