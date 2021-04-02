function compare_wall_velo(sim_dir, z_in, z_out, mu, rho, R, c_n, B_n, G_n, T, n_modes, t_steps, sol_idx)

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
w_wall_ax  = subplot(1, 2, 1); hold(w_wall_ax, 'on');
vr_wall_ax = subplot(1, 2, 2); hold(vr_wall_ax, 'on');
w_wall_hlines = []; vr_wall_hlines = [];

for ii = 1 : (t_steps + 1)
    
    % Paraview's plot-over-line data on wall-axis (y=R) ===================================
    filename = [sim_dir, '/pv_plot-over-Wallaxis_', sprintf('%06d', sol_idx(ii)), '.csv'];
    disp(['Reading ', filename]);
    
    data_interp = readmatrix(filename);
    
    x_interp  = data_interp(:, 13);
    y_interp  = data_interp(:, 14);
    z_interp  = data_interp(:, 15);
    theta_interp = atan2(y_interp, x_interp);
    
    u_interp = data_interp(:, 4);
    v_interp = data_interp(:, 5);
    w_interp = data_interp(:, 6);
    
    % Cartesian to polar transformation
    vr_interp = cos(theta_interp) .* u_interp + sin(theta_interp) .* v_interp;
    
    % Verify angular velocity is ~zero
    vt_interp = -sin(theta_interp) .* u_interp + cos(theta_interp) .* v_interp;
    
    t = (ii - 1) * dt;
    
    w_exact  = zeros(1, length(z_exact));         % axial velocity
    vr_exact = zeros(1, length(z_exact));         % radial velocity
    
    for k = 2 : n_modes
        
        n = k - 1;
        Omega_n  = R * sqrt(rho * n * omega / mu);   % Womersley number
        Lambda_n = 1j^1.5 * Omega_n;

        w_exact = w_exact + B_n(k) / ( rho * c_n(k) ) * ...
            ( 1 - G_n(k) * besselj(0, Lambda_n) / besselj(0, Lambda_n) ) * ...
            exp(1j * n * omega * (t - z_exact / c_n(k)) );
        vr_exact = vr_exact + B_n(k) * 1j * n * omega * R / ( 2 * rho * c_n(k)^2 ) .* ...
            ( 1 - G_n(k) * 2 * besselj(1, Lambda_n) / (Lambda_n * besselj(0, Lambda_n)) ) * ...
            exp(1j * n * omega * (t - z_exact / c_n(k)) );
        
    end
    
    h = plot(w_wall_ax, z_exact, real(w_exact), 'Color', colors(ii, :), 'Linestyle', '-', 'LineWidth', 1);
    w_wall_hlines = [w_wall_hlines, h];
    plot(w_wall_ax, z_interp, w_interp, 'Color', colors(ii, :), 'Linestyle', '--', 'LineWidth', 1);
    
    h = plot(vr_wall_ax, z_exact, real(vr_exact), 'Color', colors(ii, :), 'Linestyle', '-', 'LineWidth', 1);
    vr_wall_hlines = [vr_wall_hlines, h];
    plot(vr_wall_ax, z_interp, vr_interp, 'Color', colors(ii, :), 'Linestyle', '--', 'LineWidth', 1);
    
end

w_wall_ax.XLim = [z_in, z_out]; vr_wall_ax.XLim = [z_in, z_out];
w_wall_ax.YLim = [-4, 4];       vr_wall_ax.YLim = [-8e-3, 8e-3];

grid(w_wall_ax, 'minor'); grid(vr_wall_ax, 'minor');

xlabel(w_wall_ax,  '{\boldmath$z$} \bf{(cm)}', 'interpreter', 'latex', ...
       'FontName', 'Helvetica', 'FontSize', 12, 'FontWeight', 'Bold');
xlabel(vr_wall_ax, '{\boldmath$z$} \bf{(cm)}', 'interpreter', 'latex', ...
       'FontName', 'Helvetica', 'FontSize', 12, 'FontWeight', 'Bold');

ylabel(w_wall_ax,  '{\boldmath$v_z(r=R, z, t)$} \bf{(cm/s)}', 'interpreter', 'latex', ...
      'FontName', 'Helvetica', 'FontSize', 12, 'FontWeight', 'Bold');
ylabel(vr_wall_ax, '{\boldmath$v_r(r=R, z, t)$} \bf{(cm/s)}', 'interpreter', 'latex', ...
       'FontName', 'Helvetica', 'FontSize', 12, 'FontWeight', 'Bold');

set(w_wall_ax, 'Box', 'on', 'TickDir', 'out', ...
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

set(vr_wall_ax, 'Box', 'on', 'TickDir', 'out', ...
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
    
axis(w_wall_ax, 'square'); axis(vr_wall_ax, 'square');

legend(w_wall_hlines,  t_labs, 'interpreter', 'latex', 'NumColumns', t_steps + 1, 'Box', 'off', ...
       'Position', [0.4, 0.1, 0.2, 0.2], 'Units', 'normalized');
legend(vr_wall_hlines, t_labs, 'interpreter', 'latex', 'NumColumns', t_steps + 1, 'Box', 'off', ...
       'Position', [0.4, 0.1, 0.2, 0.2], 'Units', 'normalized');

set(gcf, 'WindowState', 'fullscreen');
print(gcf,  [sim_dir, '/exact-numer_wall-velo.pdf'], '-dpdf', '-r0', '-fillpage');
