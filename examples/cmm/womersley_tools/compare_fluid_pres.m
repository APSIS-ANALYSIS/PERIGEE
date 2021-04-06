function compare_fluid_pres(sim_dir, z_in, z_out, p0, mu, R, c_n, B_n, Q_n, T, n_modes, t_steps, sol_idx)

conversion = 1333.2;

colors = [     0, 0.4470, 0.7410; 0.8500, 0.3250, 0.0980; ...
          0.4940, 0.1840, 0.5560; 0.4660, 0.6740, 0.1880; ...
          0.3010, 0.7450, 0.9330; 0.6350, 0.0780, 0.1840];

linewidths = [1, 1.6];
linestyles = {'--', ':'};

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
subplot(1, 2, 1); hold on;
h_lines = [];

for ii = 1 : (t_steps + 1)
    
    % Paraview's plot-over-line data on z-axis ===================================
    num_sim = length(sim_dir);
    z_interp = cell(1, num_sim); p_interp = cell(1, num_sim);
    
    for jj = 1 : num_sim
        filename = [sim_dir{jj}, '/pv_plot-over-Zaxis_', sprintf('%06d', sol_idx{jj}(ii)), '.csv'];
        disp(['Reading ', filename]);

        data_interp = readmatrix(filename);

        z_interp{jj} = data_interp(:, 15);
        p_interp{jj} = data_interp(:,  2);
    end
    
    t = (ii - 1) * dt;
    
    % Initialize with Poiseuille solution (0th mode)
    p_exact = p0 - 8 * mu * Q_n(1) / (pi * R^4) * z_exact;     % equivalent to p0 + k0*z
    
    for k = 2 : n_modes
        
        n = k - 1;
        p_exact = p_exact + B_n(k) * exp(1j * n * omega * (t - z_exact / c_n(k)) );
        
    end
    
    h = plot(z_exact, real(p_exact) / conversion, 'Color', colors(ii, :), 'LineWidth', 1, 'Linestyle', '-');
    h_lines = [h_lines, h];
    
    for jj = 1 : num_sim
        plot(z_interp{jj}, p_interp{jj} / conversion, 'Color', colors(ii, :), ...
             'LineWidth', linewidths(jj), 'Linestyle', linestyles{jj});
    end

end

hXLabel = xlabel('{\boldmath$z$} \bf{(cm)}', 'interpreter', 'latex');
hYLabel = ylabel('{\boldmath$P(z, t)$} \bf{ (mm Hg)}', 'interpreter', 'latex');
set([hXLabel, hYLabel], 'FontName', 'Helvetica', 'FontSize', 12, 'FontWeight', 'bold');

set(gca, 'Box', 'on', 'TickDir', 'out', ...
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
    
axis square; grid minor; xlim([z_in, z_out]); ylim([-6, 6]);

lg = legend(h_lines, t_labs, 'interpreter', 'latex', 'NumColumns', 3, 'Box', 'off');
set(lg, 'Position', [0.4, 0.1, 0.2, 0.2], 'Units', 'normalized');

set(gcf, 'WindowState', 'fullscreen');

if num_sim > 1
    print(gcf, 'exact-numer_fluid-pressures.pdf', '-dpdf', '-r0', '-fillpage');
else
    print(gcf, [sim_dir{1}, '/exact-numer_fluid-pressures.pdf'], '-dpdf', '-r0', '-fillpage');
end