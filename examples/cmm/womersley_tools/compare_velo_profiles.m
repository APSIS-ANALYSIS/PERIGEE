function compare_velo_profiles(sim_dir, solver, sim_labels, z_coord, mu, rho, R, c_n, B_n, Q_n, G_n, T, n_modes, t_steps, start_step, stop_step, sol_idx)

colors = [0.918, 0.235, 0.325; 0, 0, 0.545];
linewidths = [1, 1.6];
linestyles = {'--', ':'};

dt = T / t_steps;
      
omega = 2 * pi / T;                                 % base angular frequency

x = -R : 0.01 : R;  % radius vector
r_exact = abs(x);

w_fig = figure;   v_fig = figure;
w_lim = [-8, 26]; v_lim = [-8e-3, 7e-3];

% sim_steps = start_step + (0 : t_steps) * (stop_step - start_step) / t_steps;

for ii = 1 : (t_steps + 1)
    

%     % v1: Assemble all nodal solutions on z=7.5 plane ===================================
%     filename = [sim_dir, '/SOL_9', sprintf('%08d', sim_steps(ii)), '_velo_z7d5.txt'];
%     disp(['Reading ', filename]);
%     
%     velo_numer = readmatrix(filename);
%     
%     x_numer = velo_numer(:, 1);
%     y_numer = velo_numer(:, 2);
%     theta_numer = atan2(y_numer, x_numer);
%     r_numer = sqrt(x_numer.^2 + y_numer.^2);
%     
%     % Flip the radius sign for y < 0
%     r_numer(y_numer < 0) = -r_numer(y_numer < 0);
%     
%     u_numer = velo_numer(:, 4);
%     v_numer = velo_numer(:, 5);
%     w_numer = velo_numer(:, 6);
%     
%     % Cartesian to polar transformation
%     vr_numer = cos(theta_numer) .* u_numer + sin(theta_numer) .* v_numer;
%     
%     % Verify angular velocity is ~zero
%     vt_numer = -sin(theta_numer) .* u_numer + cos(theta_numer) .* v_numer;

    % v2: Paraview's plot-over-line data on y-axis in z=7.5 plane ===================================
    num_sim = length(sim_dir); 
    
    x_interp  = cell(1, num_sim); y_interp  = cell(1, num_sim); theta_interp = cell(1, num_sim);
    u_interp  = cell(1, num_sim); v_interp  = cell(1, num_sim); w_interp = cell(1, num_sim);
    vr_interp = cell(1, num_sim); vt_interp = cell(1, num_sim);
    
    for jj = 1 : num_sim
        filename = [sim_dir{jj}, '/pv_plot-over-Yaxis_', sprintf('%06d', sol_idx{jj}(ii)), '.csv'];
        disp(['Reading ', filename]);

        velo_interp = readmatrix(filename);
        
        if strcmp(solver{jj}, 'pg')
            x_interp{jj} = velo_interp(:, 13);
            y_interp{jj} = velo_interp(:, 14);
            
            u_interp{jj} = velo_interp(:, 4);
            v_interp{jj} = velo_interp(:, 5);
            w_interp{jj} = velo_interp(:, 6);
        
        elseif strcmp(solver{jj}, 'sv')
            x_interp{jj} = velo_interp(:, 26);
            y_interp{jj} = velo_interp(:, 27);
            
            u_interp{jj} = velo_interp(:, 17);
            v_interp{jj} = velo_interp(:, 18);
            w_interp{jj} = velo_interp(:, 19);
            
        else
            disp('Unknown solver');
        end

        theta_interp{jj} = atan2(y_interp{jj}, x_interp{jj});
        
        % Cartesian to polar transformation
        vr_interp{jj} =  cos(theta_interp{jj}) .* u_interp{jj} + sin(theta_interp{jj}) .* v_interp{jj};

        % Verify angular velocity is ~zero
        vt_interp{jj} = -sin(theta_interp{jj}) .* u_interp{jj} + cos(theta_interp{jj}) .* v_interp{jj};
    end
    
    t = (ii - 1) * dt;
    w_ax = subplot((t_steps + 1) / 2, 2, ii, 'Parent', w_fig); hold(w_ax, 'on');
    v_ax = subplot((t_steps + 1) / 2, 2, ii, 'Parent', v_fig); hold(v_ax, 'on');
    
    ax_all = {w_ax, v_ax};
    
    w_exact = 2 * Q_n(1) * (R^2 - r_exact.^2) / (pi * R^4);         % axial velocity

    vr_exact = zeros(1, length(r_exact));                           % radial velocity

    for k = 2 : n_modes
        
        n = k - 1;
        Omega_n  = R * sqrt(rho * n * omega / mu);   % Womersley number
        Lambda_n = 1j^1.5 * Omega_n;

        w_exact = w_exact + B_n(k) / ( rho * c_n(k) ) * ...
            ( 1 - G_n(k) * besselj(0, Lambda_n * r_exact / R) / besselj(0, Lambda_n) ) * ...
            exp(1j * n * omega * (t - z_coord / c_n(k)) );
        vr_exact = vr_exact + B_n(k) * 1j * n * omega * R / ( 2 * rho * c_n(k)^2 ) .* ...
            ( r_exact / R - G_n(k) * 2 * besselj(1, Lambda_n * r_exact / R) / (Lambda_n * besselj(0, Lambda_n)) ) * ...
            exp(1j * n * omega * (t - z_coord / c_n(k)) );
        
    end
    
    plot(w_ax, [0, 0], [-R, R], 'Color', [0.75, 0.75, 0.75] );               % grey
    plot(w_ax, real(w_exact), x, 'Color', colors(1, :), 'LineWidth', 1, 'Linestyle', '-');
    % plot(w_ax, w_numer, r_numer, 'Color', colors(2, :), 'Linestyle', 'None', 'Marker', 'o', 'MarkerSize', 3);
    
    for jj = 1 : num_sim
        plot(w_ax, w_interp{jj}, y_interp{jj}, 'Color', colors(2, :), ...
             'LineWidth', linewidths(jj), 'Linestyle', linestyles{jj});
    end
    
    set(w_ax, 'Box', 'on', 'TickDir', 'out', ...
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
    
    plot(v_ax, [0, 0], [-R, R], 'Color', [0.75, 0.75, 0.75] );               % grey
    plot(v_ax, real(vr_exact), x, 'Color', colors(1, :), 'LineWidth', 1, 'Linestyle', '-');
    % plot(v_ax, vr_numer, r_numer, 'Color', colors(2, :), 'Linestyle', 'None', 'Marker', 'o', 'MarkerSize', 3);
    
    for jj = 1 : num_sim
        plot(v_ax, vr_interp{jj}, y_interp{jj}, 'Color', colors(2, :), ...
             'LineWidth', linewidths(jj), 'Linestyle', linestyles{jj});
    end
    
    set(v_ax, 'Box', 'on', 'TickDir', 'out', ...
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
    
    w_ax.XLim  = w_lim;  v_ax.XLim  = v_lim;
    w_ax.YLim = [-R, R+0.001]; v_ax.YLim = [-R, R+0.001];
          
    ylabel(w_ax, '{\boldmath$y$} \bf{(cm)}', 'interpreter', 'latex', ...
           'FontName', 'Helvetica', 'FontSize', 12, 'FontWeight', 'Bold');
    ylabel(v_ax, '{\boldmath$y$} \bf{(cm)}', 'interpreter', 'latex', ...
           'FontName', 'Helvetica', 'FontSize', 12, 'FontWeight', 'Bold');
    
    axis(w_ax, 'square'); axis(v_ax, 'square');
    grid(w_ax, 'minor');  grid(v_ax, 'minor');
    
    if ii == 1
        w_lg = legend(w_ax, ['Analytical', sim_labels], 'NumColumns', num_sim + 1, 'Box', 'off');
        v_lg = legend(v_ax, ['Analytical', sim_labels], 'NumColumns', num_sim + 1, 'Box', 'off');
    end
    
    for jj = 1 : length(ax_all)
        if ii == 1
            text(w_ax, 0.6, 0.9, '{\boldmath$t$ = $0$}', ...
                 'interpreter', 'latex', 'Units', 'normalized', ...
                 'FontName', 'Helvetica', 'FontSize', 12, 'FontWeight', 'Bold');
            text(v_ax, 0.1, 0.9, '{\boldmath$t$ = $0$}', ...
                 'interpreter', 'latex', 'Units', 'normalized', ...
                 'FontName', 'Helvetica', 'FontSize', 12, 'FontWeight', 'Bold');
        elseif ii == 2
            text(w_ax, 0.6, 0.9, ['{\boldmath$t$ = $T$/}\bf{', num2str(t_steps), '}'], ...
                 'interpreter', 'latex', 'Units', 'normalized', ...
                 'FontName', 'Helvetica', 'FontSize', 12, 'FontWeight', 'Bold');
            text(v_ax, 0.1, 0.9, ['{\boldmath$t$ = $T$/}\bf{', num2str(t_steps), '}'], ...
                 'interpreter', 'latex', 'Units', 'normalized', ...
                 'FontName', 'Helvetica', 'FontSize', 12, 'FontWeight', 'Bold');
        elseif ii == t_steps + 1
            text(w_ax, 0.6, 0.9, '{\boldmath$t$ = $T$}', ...
                 'interpreter', 'latex', 'Units', 'normalized', ...
                 'FontName', 'Helvetica', 'FontSize', 12, 'FontWeight', 'Bold');
            text(v_ax, 0.1, 0.9, '{\boldmath$t$ = $T$}', ...
                 'interpreter', 'latex', 'Units', 'normalized', ...
                 'FontName', 'Helvetica', 'FontSize', 12, 'FontWeight', 'Bold');
        else
            text(ax_all{jj}, 0.6, 0.9, ['{\boldmath$t$ = }\bf{', num2str(ii-1), '}{\boldmath$T$/}', num2str(t_steps)], ...
                 'interpreter', 'latex', 'Units', 'normalized', ...
                 'FontName', 'Helvetica', 'FontSize', 12, 'FontWeight', 'Bold');
        end
    end
end

set(w_lg, 'Position', [0.4, -0.08, 0.2, 0.2], 'Units', 'normalized');
set(v_lg, 'Position', [0.4, -0.08, 0.2, 0.2], 'Units', 'normalized');

sgtitle(w_fig, '{\boldmath$v_z(r, z=L/2, t)$} \bf{(cm/s)}', ...
        'interpreter', 'latex', 'FontName', 'Helvetica', 'FontSize', 12, 'FontWeight', 'bold');
sgtitle(v_fig, '{\boldmath$v_r(r, z=L/2, t)$} \bf{(cm/s)}', ...
        'interpreter', 'latex', 'FontName', 'Helvetica', 'FontSize', 12, 'FontWeight', 'bold');

set(w_fig, 'WindowState','fullscreen');
set(v_fig, 'WindowState','fullscreen');

if num_sim > 1
    print(v_fig, 'exact-numer_radial-velo-profiles.pdf', '-dpdf', '-r0', '-fillpage');
    print(w_fig, 'exact-numer_axial-velo-profiles.pdf',  '-dpdf', '-r0', '-fillpage');
else
    print(v_fig, [sim_dir{1}, '/exact-numer_radial-velo-profiles.pdf'], '-dpdf', '-r0', '-fillpage');
    print(w_fig, [sim_dir{1}, '/exact-numer_axial-velo-profiles.pdf'],  '-dpdf', '-r0', '-fillpage');
end