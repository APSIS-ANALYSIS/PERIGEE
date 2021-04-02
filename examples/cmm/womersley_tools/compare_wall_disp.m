function compare_wall_disp(sim_dir, z_in, z_out, rho, R, c_n, g_n, B_n, G_n, T, n_modes, t_steps, sol_idx)

colors = [     0, 0.4470, 0.7410; 0.8500, 0.3250, 0.0980; ...
          0.4940, 0.1840, 0.5560; 0.4660, 0.6740, 0.1880; ...
          0.3010, 0.7450, 0.9330; 0.6350, 0.0780, 0.1840];

dt = T / t_steps;
      
omega = 2 * pi / T;                                 % base angular frequency

t_labs = cell(1, t_steps + 1);
t_labs{1} = '{\boldmath$t$ = $0$}';
t_labs{2} = ['{\boldmath$t$ = $T$/}\bf{', num2str(t_steps), '}'];
for ii = 3 : t_steps
    t_labs{ii} = ['{\boldmath$t$ = }\bf{', num2str(ii-1), '}{\boldmath$T$/}\bf{', num2str(t_steps), '}'];
end
t_labs{t_steps + 1} = '{\boldmath$t$ = $T$}';

z_exact = linspace(z_in, z_out, 100);

figure;
xi_ax  = subplot(1, 2, 1); hold(xi_ax, 'on');
eta_ax = subplot(1, 2, 2); hold(eta_ax, 'on');
xi_hlines = []; eta_hlines = [];

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

xi_ax.XLim = [z_in, z_out]; eta_ax.XLim = [z_in, z_out];
xi_ax.YLim = [-0.7, 0.7];   eta_ax.YLim = [-1.5e-3, 1.5e-3];

xlabel(xi_ax,  '{\boldmath$z$} \bf{(cm)}', 'interpreter', 'latex', ...
       'FontName', 'Helvetica', 'FontSize', 12, 'FontWeight', 'Bold');
xlabel(eta_ax, '{\boldmath$z$} \bf{(cm)}', 'interpreter', 'latex', ...
       'FontName', 'Helvetica', 'FontSize', 12, 'FontWeight', 'Bold');

ylabel(xi_ax,  '{\boldmath$u_z(z, t)$} \bf{(cm)}', 'interpreter', 'latex', ...
      'FontName', 'Helvetica', 'FontSize', 12, 'FontWeight', 'Bold');
ylabel(eta_ax, '{\boldmath$u_r(z, t)$} \bf{(cm)}', 'interpreter', 'latex', ...
       'FontName', 'Helvetica', 'FontSize', 12, 'FontWeight', 'Bold');

set(xi_ax, 'Box', 'on', 'TickDir', 'out', ...
        'TickLength'  , [.02 .02], ...
        'XMinorTick'  , 'on' , ...
        'YMinorTick'  , 'on' , ...
        'YGrid'       , 'on' , ...
        'XGrid'       , 'on' , ...
        'XColor'      , [0 0 0 ], ...
        'YColor'      , [0 0 0 ], ...
        'LineWidth'   , 1, ...
        'FontSize', 12, ...
        'FontWeight', 'Bold');

set(eta_ax, 'Box', 'on', 'TickDir', 'out', ...
        'TickLength'  , [.02 .02], ...
        'XMinorTick'  , 'on' , ...
        'YMinorTick'  , 'on' , ...
        'YGrid'       , 'on' , ...
        'XGrid'       , 'on' , ...
        'XColor'      , [0 0 0 ], ...
        'YColor'      , [0 0 0 ], ...
        'LineWidth'   , 1, ...
        'FontSize', 12, ...
        'FontWeight', 'Bold');
    
axis(xi_ax, 'square'); axis(eta_ax, 'square');
grid(xi_ax, 'minor');  grid(eta_ax, 'minor');

legend(xi_hlines,  t_labs, 'interpreter', 'latex', 'NumColumns', t_steps + 1, 'Box', 'off', ...
       'Position', [0.4, 0.1, 0.2, 0.2], 'Units', 'normalized');
legend(eta_hlines, t_labs, 'interpreter', 'latex', 'NumColumns', t_steps + 1, 'Box', 'off', ...
       'Position', [0.4, 0.1, 0.2, 0.2], 'Units', 'normalized');

set(gcf, 'WindowState', 'fullscreen');
print(gcf, [sim_dir, '/exact-numer_wall-disp.pdf'], '-dpdf', '-r0', '-fillpage');
