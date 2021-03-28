function compare_wall_disp(sim_dir, z_in, z_out, rho, R, c_n, g_n, B_n, G_n, T, n_modes, t_steps, sol_idx)

colors = [     0, 0.4470, 0.7410; 0.8500, 0.3250, 0.0980; ...
          0.4940, 0.1840, 0.5560; 0.4660, 0.6740, 0.1880; ...
          0.3010, 0.7450, 0.9330; 0.6350, 0.0780, 0.1840];

dt = T / t_steps;
      
omega = 2 * pi / T;                                 % base angular frequency

t_labs = cell(1, t_steps + 1);
t_labs{1} = '$t$ = $0$';
t_labs{2} = ['$t$ = $T$ / ', num2str(t_steps)];
for ii = 3 : t_steps
    t_labs{ii} = ['$t$ = ', num2str(ii-1), ' $T$ / ', num2str(t_steps)];
end
t_labs{t_steps + 1} = '$t$ = $T$';

z_exact = linspace(z_in, z_out, 100);

xi_fig  = figure;  xi_ax = axes;  hold on; xi_hlines  = [];
eta_fig = figure; eta_ax = axes;  hold on; eta_hlines = [];


for ii = 1 : (t_steps + 1)
    
    % Paraview's plot-over-line data on wall-axis (y=R) ===================================
    filename = [sim_dir, '/pv_plot-over-Wallaxis_', sprintf('%06d', sol_idx(ii)), '.csv'];
    disp(['Reading ', filename]);
    
    data_interp = readmatrix(filename);
    
    x_interp  = data_interp(:, 13);
    y_interp  = data_interp(:, 14);
    z_interp  = data_interp(:, 15);
    theta_interp = atan2(y_interp, x_interp);
    
    ux_interp = data_interp(:, 7);          % x-displacement
    uy_interp = data_interp(:, 8);          % y-displacement
    xi_interp = data_interp(:, 9);          % z-displacement
       
    % Cartesian to polar transformation
    eta_interp = cos(theta_interp) .* ux_interp + sin(theta_interp) .* uy_interp;
%     r0_numer = sqrt( x_interp.^2 + y_interp.^2 );
%     r_numer  = sqrt( (x_interp + ux_interp).^2 + (y_interp + uy_interp).^2 );
%     eta_interp = r_numer - r0_numer;
    
    % Verify angular displacement is ~zero
    phi_interp = -sin(theta_interp) .* ux_interp + cos(theta_interp) .* uy_interp;
    
    t = (ii - 1) * dt;
    
    xi_exact  = zeros(1, length(z_exact));        % axial displacement
    eta_exact = zeros(1, length(z_exact));        % radial displacement
    
    for k = 2 : n_modes
        
        n = k - 1;
        xi_exact  = xi_exact + 1j * B_n(k) / ( rho * c_n(k) * n * omega ) * ( G_n(k) - 1 ) * ...
                    exp(1j * n * omega * (t - z_exact / c_n(k)) );
        eta_exact = eta_exact + B_n(k) * R / ( 2 * rho * c_n(k)^2 ) * ( 1 - G_n(k) * g_n(k) ) * ...
                    exp(1j * n * omega * (t - z_exact / c_n(k)) ); 
        
    end
    
    h = plot(xi_ax, z_exact, real(xi_exact), 'Color', colors(ii, :), 'Linestyle', '-', 'LineWidth', 1);
    xi_hlines = [xi_hlines, h];
    plot(xi_ax, z_interp, xi_interp, 'Color', colors(ii, :), 'Linestyle', '--', 'LineWidth', 1);
    
    h = plot(eta_ax, z_exact, real(eta_exact), 'Color', colors(ii, :), 'Linestyle', '-', 'LineWidth', 1);
    eta_hlines = [eta_hlines, h];
    plot(eta_ax, z_interp, eta_interp, 'Color', colors(ii, :), 'Linestyle', '--', 'LineWidth', 1);
    
end

grid(xi_ax, 'minor'); grid(eta_ax, 'minor');

xlabel(xi_ax,  'z (cm)', 'interpreter', 'latex');
xlabel(eta_ax, 'z (cm)', 'interpreter', 'latex');

ylabel(xi_ax,  'cm', 'interpreter', 'latex');
ylabel(eta_ax, 'cm', 'interpreter', 'latex');

title(xi_ax,  'Axial Wall Disp $\xi(z, t)$',   'interpreter', 'latex');
title(eta_ax, 'Radial Wall Disp $\eta(z, t)$', 'interpreter', 'latex');

legend(xi_hlines,  t_labs, 'interpreter', 'latex', 'NumColumns', 3, 'Location', 'best');
legend(eta_hlines, t_labs, 'interpreter', 'latex', 'NumColumns', 3, 'Location', 'best');

saveas(xi_fig,  [sim_dir, '/exact-numer_axial-wall-disp.png']);
saveas(eta_fig, [sim_dir, '/exact-numer_radial-wall-disp.png']);