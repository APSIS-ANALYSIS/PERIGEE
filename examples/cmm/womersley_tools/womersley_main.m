% Deformable Womersley solution
% Performs FFT on an inflow profile
% Q(z, t) = B_n * pi * R^2 / (rho * c) * (1 - G*g) * exp( inw * (t - z/c) )
% Q(0, t) = B_n * pi * R^2 / (rho * c) * (1 - G*g) * exp( inwt )
% and computes the following quantities
% pressure:                 p(r, z, t)
% axial velocity:           w(r, z, t); radial velocity:          v(r, z, t)
% axial  wall displacement: xi(z, t);   radial wall displacement: eta(z, t)

close all; clear; clc;

T = 1.1;                            % period (s)

n_modes = 2;                        % num Fourier modes (including steady 0th mode)

% Fluid properties
mu = 0.04;                          % dynamic viscosity (dyn * s / cm2)
rho = 1.0;                          % density (g / cm3)
nu = mu / rho;                      % kinematic viscosity (stokes)
p0 = 0.0; % 133322;                 % mean inlet pressure (dyn / cm2)


% Wall properties
R     = 0.3;                                                    % tube radius (cm)
rho_s = 1.0;                                                    % wall density (g / cm3)
nu_s  = 0.5;                                                    % wall Poisson's ratio
h_s   = 0.06; % 0.03;                                           % wall thickness (cm)
E_s   = 1.0e4 * 13.3^2 / (h_s * (0.2 * R)^(-0.4)); % 9863400;   % wall Young's modulus (dyn/cm2) 

[c_n, gamma_n, g_n] = wave_speed(mu, rho, T, R, n_modes, rho_s, nu_s, h_s, E_s);

% Flow rate (mL/s) measured in a pig's main pulmonary artery
flow = [ 33.42, 56.19, 73.697, 96.721, 139.85, 164.46, 177.44, 196.25, ...
         198.77, 184.72, 162.09, 131.85, 91.057, 75.404, 62.991, 32.539, ...
         21.865, 28.182, 23.896, 19.457, 19.911, 13.432, 5.284, -1.0584 ];

[B_n, Q_n, G_n] = compute_B(flow / 50, rho, nu_s, c_n, gamma_n, g_n, T, R, n_modes);

% Plot velocity profiles over time down the length of the tube
z_steps = 5;        % num intervals in [0, first wavelength]
t_steps = 5;        % num intervals in [0, T]
plot_fluid_velo(mu, rho, R, c_n, B_n, Q_n, G_n, T, n_modes, z_steps, t_steps);

% Plot wall displacement over time down the length of the tube
plot_wall_disp(rho, R, c_n, g_n, B_n, G_n, T, n_modes, t_steps);

% Plot fluid pressure over time down the length of the tube
plot_fluid_pres(p0, mu, R, c_n, B_n, Q_n, T, n_modes, t_steps);