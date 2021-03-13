% Compute the wave speed c as a function of the wall properties and 
% Womersley number

function [c_n, gamma_n, g_n] = wave_speed(mu, rho, T, R, n_modes, rho_s, nu_s, h_s, E_s)

omega = 2 * pi / T;                             % base angular frequency
assignin('base', 'omega', omega);

c_n = zeros(1, n_modes);                     % (complex) viscous wave speed
c_R = zeros(1, n_modes);                     % real wave speed (dispersion coefficient)
c_I = zeros(1, n_modes);                     % imag wave speed (attenuation coefficient)

g_n = zeros(1, n_modes);
gamma_n = zeros(1, n_modes);                   % solution to frequency eqn
 
for k = 2 : n_modes
    n = k - 1;
    Omega_n  = R * sqrt(rho * n * omega / mu);   % Womersley number
    Lambda_n = 1j^1.5 * Omega_n;
    g_n(k) = 2 * besselj(1, Lambda_n) / ( Lambda_n * besselj(0, Lambda_n) );
    
    % Define coefficients of the quadratic frequency equation in gamma as
    % p * gamma^2 + q * gamma + s = 0
    p = (g_n(k) - 1) * (nu_s^2 - 1);
    q = rho_s * h_s / (rho * R) * (g_n(k) - 1) + (2 * nu_s - 0.5) * g_n(k) - 2;
    s = 2 * rho_s * h_s / (rho * R) + g_n(k);
    
    gamma1 = -q + sqrt(q^2 - 4 * p * s) / (2 * p);
    gamma2 = -q - sqrt(q^2 - 4 * p * s) / (2 * p);
    
    gamma_n(k) = max(gamma1, gamma2);
    
    % Inviscid wave speed (Moens-Korteweg)
    c0 = sqrt(E_s * h_s / (2 * rho * R) );
    
    c_n(k)   = c0 * sqrt(2 / ( (1 - nu_s^2) * gamma_n(k) ) );
    c_R(k) = 1 / real(1 / c_n(k));
    c_I(k) = 1 / imag(1 / c_n(k));

end

assignin('base', 'c_R', c_R);
assignin('base', 'c_I', c_I);