function compare_flow_pres(sim_dir, z_in, z_out, inlet_data, outlet_data, p0, mu, rho, R, c_n, B_n, Q_n, G_n, g_n, T, n_modes, sol_idx)

conversion = 1333.2;

colors = [0.918, 0.235, 0.325; 0, 0, 0.545];
linewidths = [1, 1.6];
linestyles = {'--', ':'};
      
omega = 2 * pi / T;                                 % base angular frequency

t_exact = linspace(0, T, 100)';

% Analytical solution
q_in_exact  = Q_n(1) * ones(length(t_exact), 1);
q_out_exact = Q_n(1) * ones(length(t_exact), 1);
p_in_exact  = ( p0 - 8 * mu * Q_n(1) / (pi * R^4) * z_in  ) * ones(length(t_exact), 1);
p_out_exact = ( p0 - 8 * mu * Q_n(1) / (pi * R^4) * z_out ) * ones(length(t_exact), 1);

for k = 2 : n_modes
    
    n = k - 1;
    q_in_exact  = q_in_exact + B_n(k) * pi * R^2 / (rho * c_n(k)) * (1 - G_n(k) * g_n(k)) * ...
                  exp(1j * n * omega * (t_exact - z_in / c_n(k)) );
    q_out_exact = q_out_exact + B_n(k) * pi * R^2 / (rho * c_n(k)) * (1 - G_n(k) * g_n(k)) .* ...
                  exp(1j * n * omega * (t_exact - z_out / c_n(k)) );
    p_in_exact  = p_in_exact  + B_n(k) * exp(1j * n * omega * (t_exact -  z_in / c_n(k)) );
    p_out_exact = p_out_exact + B_n(k) * exp(1j * n * omega * (t_exact - z_out / c_n(k)) );
    
end

num_sim = length(sim_dir); 
t_numer = cell(1, num_sim); 
q_in_numer = cell(1, num_sim); q_out_numer = cell(1, num_sim);
p_in_numer = cell(1, num_sim); p_out_numer = cell(1, num_sim);

for ii = 1 : num_sim
    t_numer{ii} = inlet_data{ii}(sol_idx{ii}, 2) - inlet_data{ii}(sol_idx{ii}(1), 2);
    
    q_in_numer{ii}  = -inlet_data{ii}(sol_idx{ii}, 4); % flip sign for unit normal
    q_out_numer{ii} = outlet_data{ii}(sol_idx{ii}, 4);
    
    p_in_numer{ii}  =  inlet_data{ii}(sol_idx{ii}, 5);
    p_out_numer{ii} = outlet_data{ii}(sol_idx{ii}, 5);

end

% ================== Plot analytical vs. numerical flows ==================
figure; 
subplot(1, 2, 1); hold on;
plot(t_exact, real(q_in_exact),  'Color', colors(1, :), ...
     'LineWidth', 1, 'Linestyle', '-');
plot(t_exact, real(q_out_exact), 'Color', colors(2, :), ...
     'LineWidth', 1, 'Linestyle', '-');
 
for ii = 1 : num_sim
    plot(t_numer{ii},  q_in_numer{ii}, 'Color', colors(1, :), ...
         'LineWidth', linewidths(ii), 'Linestyle', linestyles{ii});
    plot(t_numer{ii}, q_out_numer{ii}, 'Color', colors(2, :), ...
         'LineWidth', linewidths(ii), 'Linestyle', linestyles{ii});
end

hXLabel = xlabel('{\boldmath$t$} \bf{(s)}', 'Interpreter', 'Latex');
hYLabel = ylabel('{\boldmath$Q(z, t)$} \bf{ (mL/s)}', 'Interpreter', 'Latex');
set([hXLabel, hYLabel], 'FontName', 'Helvetica', 'FontSize', 12, 'FontWeight', 'bold');

set( gca, 'Box', 'on', 'TickDir'     , 'out', ...
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
    
axis square; grid minor; xlim([0, T]);

% ================  Plot analytical vs. numerical pressures ===============
subplot(1, 2, 2); hold on;
plot(t_exact, real(p_in_exact)  / conversion, 'Color', colors(1, :), ...
     'LineWidth', 1, 'Linestyle', '-');
plot(t_exact, real(p_out_exact) / conversion, 'Color', colors(2, :), ...
     'LineWidth', 1, 'Linestyle', '-');
 
for ii = 1 : num_sim
    plot(t_numer{ii},  p_in_numer{ii} / conversion, 'Color', colors(1, :), ...
         'LineWidth', linewidths(ii), 'Linestyle', linestyles{ii});
    plot(t_numer{ii}, p_out_numer{ii} / conversion, 'Color', colors(2, :), ...
         'LineWidth', linewidths(ii), 'Linestyle', linestyles{ii});
end

hXLabel = xlabel('{\boldmath$t$} \bf{(s)}', 'Interpreter', 'Latex');
hYLabel = ylabel('{\boldmath$P(z, t)$} \bf{ (mm Hg)}', 'Interpreter', 'Latex');
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
    
axis square; grid minor; xlim([0, T]); ylim([-6, 6]);
set(gca, 'FontSize', 12, 'fontWeight', 'bold');

lg = legend('Inlet / Analytical', 'Outlet / Analytical', 'Inlet / Numerical', 'Outlet / Numerical', ...
           'NumColumns', 4, 'Box', 'off');
set(lg, 'Position', [0.4, 0.1, 0.2, 0.2], 'Units', 'normalized');

set(gcf, 'WindowState','fullscreen');

if num_sim > 1
    print(gcf, 'exact-numer_cap-flows-pressures.pdf', '-dpdf', '-r0', '-fillpage');
else
    print(gcf, [sim_dir{1}, '/exact-numer_cap-flows-pressures.pdf'], '-dpdf', '-r0', '-fillpage');
end


