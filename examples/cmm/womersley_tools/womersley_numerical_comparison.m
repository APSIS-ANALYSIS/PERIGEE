% Compare analytical and numerical Womersley solutions

close all; clear; clc;

sim_dir = {'/home/ingridxlan/Documents/Ingrid/Solvers/CMM_testing/womersley_cylinder/Simulations/R0p3_L15_deformable/P1_axial_ms_5p000e-2_sv', ...
           '/home/ingridxlan/Documents/Ingrid/Solvers/CMM_testing/womersley_cylinder/Simulations/R0p3_L15_deformable/P1_axial_ms_5p000e-2', ...
           '/home/ingridxlan/Documents/Ingrid/Solvers/CMM_testing/womersley_cylinder/Simulations/R0p3_L15_deformable/P2_axial_ms_5p000e-2'};
num_sim = length(sim_dir);
solver = {'sv', 'pg', 'pg'};
sim_labels = {'svSolver', 'P1', 'P2'};      % {'Numerical'};

zoomed_views = true;

% Simulation steps
start_step = [4400, 4400, 8800];
incr_step = [40, 40, 80];
stop_step = [6600, 6600, 13200];       
num_cyc = [3, 3, 3];                % num cardiac cycles simulated

z_in  = 0;                          % z-coord of inlet face
z_out = 15;                         % z-coord of outlet face
z_half = z_out / 2;                 % z-coord halfway down the tube

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

% Compare inlet & outlet flows & pressures
inlet_data = cell(1, num_sim); outlet_data = cell(1, num_sim);
sol_idx = cell(1, num_sim);
for ii = 1 : num_sim
    
    if strcmp(solver{ii}, 'pg')
        inlet_data{ii}  = readmatrix([sim_dir{ii}, '/Outlet_000_data.txt']);
        outlet_data{ii} = readmatrix([sim_dir{ii}, '/Outlet_001_data.txt']);
    
    elseif strcmp(solver{ii}, 'sv')
        cap_pres = readmatrix([sim_dir{ii}, '/cap_pressures.txt']);
        cap_flow = readmatrix([sim_dir{ii}, '/cap_flows.txt']);
        inlet_data{ii}  = [cap_flow(:, 1 : 2),  cap_pres(:, 2)];
        outlet_data{ii} = [cap_flow(:, [1, 3]), cap_pres(:, 3)];
    
    else
        disp('Unknown solver');
    end
    
    sol_idx{ii} = (start_step(ii) : stop_step(ii)) + 1;
end

compare_flow_pres(sim_dir, solver, z_in, z_out, inlet_data, outlet_data, ...
                  p0, mu, rho, R, c_n, B_n, Q_n, G_n, g_n, T, n_modes, sol_idx, zoomed_views);

t_steps = 5;                                % number of intervals for comparison

% % Compare velocity profiles halfway down the tube
% sol_idx_all = cell(1, num_sim);
% for ii = 1 : num_sim
%     sol_idx_all{ii} = ( (num_cyc(ii) - 1) / num_cyc(ii) * stop_step(ii) / incr_step(ii)) : ...
%                       ( stop_step(ii) / incr_step(ii) );
%     sol_idx{ii} = sol_idx_all{ii}(1 : (length(sol_idx_all{ii}) - 1) / t_steps :  end);
% end
% 
% compare_velo_profiles(sim_dir, solver, sim_labels, z_half, mu, rho, R, c_n, B_n, ...
%                       Q_n, G_n, T, n_modes, t_steps, start_step, stop_step, sol_idx, zoomed_views);

% % Compare pressures down the tube
% compare_fluid_pres(sim_dir, solver, z_in, z_out, p0, mu, R, c_n, B_n, Q_n, T, ...
%                    n_modes, t_steps, sol_idx, zoomed_views);

% % Compare wall displacements down the tube
% compare_wall_disp(sim_dir, solver, z_in, z_out, rho, R, c_n, g_n, B_n, G_n, T, ...
%                   n_modes, t_steps, sol_idx, zoomed_views);
%  
% % Compare wall velocities down the tube
% compare_wall_velo(sim_dir, solver, z_in, z_out, mu, rho, R, c_n, B_n, G_n, T, ...
%                   n_modes, t_steps, sol_idx, zoomed_views);
