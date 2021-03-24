function compare_fluid_pres(sim_dir, z_in, z_out, p0, mu, R, c_n, B_n, Q_n, T, n_modes, t_steps, start_step, stop_step, sol_idx)

conversion = 1333.2;

colors = [     0, 0.4470, 0.7410; 0.8500, 0.3250, 0.0980; ...
          0.4940, 0.1840, 0.5560; 0.4660, 0.6740, 0.1880; ...
          0.3010, 0.7450, 0.9330; 0.6350, 0.0780, 0.1840];

dt = T / t_steps;
      
omega = 2 * pi / T;                                 % base angular frequency

sim_steps = start_step + (0 : t_steps) * (stop_step - start_step) / t_steps;

t_labs = cell(1, t_steps + 1);
t_labs{1} = '$t$ = $0$';
t_labs{2} = ['$t$ = $T$ / ', num2str(t_steps)];
for ii = 3 : t_steps
    t_labs{ii} = ['$t$ = ', num2str(ii-1), ' $T$ / ', num2str(t_steps)];
end
t_labs{t_steps + 1} = '$t$ = $T$';

z_exact = linspace(z_in, z_out, 100);

figure; hold on;
h_lines = [];

for ii = 1 : (t_steps + 1)
    
    % Paraview's plot-over-line data on z-axis ===================================
    filename = [sim_dir, '/pv_plot-over-Zaxis_', sprintf('%06d', sol_idx(ii)), '.csv'];
    disp(['Reading ', filename]);
    
    data_interp = readmatrix(filename);
    
    z_interp = data_interp(:, 15);
    p_interp = data_interp(:, 2);
    
    t = (ii - 1) * dt;
    
    % Initialize with Poiseuille solution (0th mode)
    p_exact = p0 - 8 * mu * Q_n(1) / (pi * R^4) * z_exact;     % equivalent to p0 + k0*z
    
    for k = 2 : n_modes
        
        n = k - 1;
        p_exact = p_exact + B_n(k) * exp(1j * n * omega * (t - z_exact / c_n(k)) );
        
    end
    
    h = plot(z_exact, real(p_exact) / conversion, 'Color', colors(ii, :), 'Linestyle', '-', 'LineWidth', 1);
    h_lines = [h_lines, h];
    plot(z_interp, p_interp / conversion, 'Color', colors(ii, :), 'Linestyle', '--', 'LineWidth', 1);

end

xlim([z_in, z_out]); grid('minor'); 
xlabel('z (cm)', 'interpreter', 'latex');
ylabel('mm Hg', 'interpreter', 'latex')
title('Fluid Pressure $P(z, t)$', 'interpreter', 'latex');
legend(h_lines, t_labs, 'interpreter', 'latex', 'NumColumns', 3, 'Location', 'best');

saveas(gcf, [sim_dir, '/exact-numer_fluid-pressures.png'])