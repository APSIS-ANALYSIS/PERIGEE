function compare_velo_profiles(sim_dir, z_coord, mu, rho, R, c_n, B_n, Q_n, G_n, T, n_modes, t_steps, start_step, stop_step, sol_idx)

colors = [     0, 0.4470, 0.7410; 0.8500, 0.3250, 0.0980; ...
          0.4940, 0.1840, 0.5560; 0.4660, 0.6740, 0.1880; ...
          0.3010, 0.7450, 0.9330; 0.6350, 0.0780, 0.1840];

dt = T / t_steps;
      
omega = 2 * pi / T;                                 % base angular frequency

x = -R : 0.01 : R;  % radius vector
r_exact = abs(x);

w_fig = figure;   v_fig = figure;
w_lim = [-5, 25]; v_lim = [-7e-3, 7e-3];

sim_steps = start_step + (0 : t_steps) * (stop_step - start_step) / t_steps;

for ii = 1 : (t_steps + 1)
    

    % v1: Assemble all nodal solutions on z=7.5 plane ===================================
    filename = [sim_dir, '/SOL_9', sprintf('%08d', sim_steps(ii)), '_velo_z7d5.txt'];
    disp(['Reading ', filename]);
    
    velo_numer = readmatrix(filename);
    
    x_numer = velo_numer(:, 1);
    y_numer = velo_numer(:, 2);
    theta_numer = atan2(y_numer, x_numer);
    r_numer = sqrt(x_numer.^2 + y_numer.^2);
    
    % Flip the radius sign for y < 0
    r_numer(y_numer < 0) = -r_numer(y_numer < 0);
    
    u_numer = velo_numer(:, 4);
    v_numer = velo_numer(:, 5);
    w_numer = velo_numer(:, 6);
    
    % Cartesian to polar transformation
    vr_numer = cos(theta_numer) .* u_numer + sin(theta_numer) .* v_numer;
    
    % Verify angular velocity is ~zero
    vt_numer = -sin(theta_numer) .* u_numer + cos(theta_numer) .* v_numer;

    % v2: Paraview's plot-over-line data on y-axis in z=7.5 plane ===================================
    filename = [sim_dir, '/pv_plot-over-line_', sprintf('%06d', sol_idx(ii)), '.csv'];
    disp(['Reading ', filename]);
    
    velo_interp = readmatrix(filename);
    
    x_interp = velo_interp(:, 13);
    y_interp = velo_interp(:, 14);
    theta_interp = atan2(y_interp, x_interp);
    
    u_interp = velo_interp(:, 4);
    v_interp = velo_interp(:, 5);
    w_interp = velo_interp(:, 6);

    % Cartesian to polar transformation
    vr_interp = cos(theta_interp) .* u_interp + sin(theta_interp) .* v_interp;
    
    % Verify angular velocity is ~zero
    vt_interp = -sin(theta_interp) .* u_interp + cos(theta_interp) .* v_interp;
    
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
    
    plot(w_ax, real(w_exact), x / R, 'Color', colors(1, :), 'Linestyle', '-', 'LineWidth', 1);
    % plot(w_ax, w_numer, r_numer / R, 'Color', colors(2, :), 'Linestyle', 'None', 'Marker', 'o', 'MarkerSize', 3);
    plot(w_ax, w_interp, y_interp / R, 'Color', colors(2, :), 'Linestyle', '--', 'LineWidth', 1);
    plot(w_ax, [0, 0], [-1, 1], 'Color', [0.75, 0.75, 0.75] );               % grey

    plot(v_ax, real(vr_exact), x / R, 'Color', colors(1, :), 'Linestyle', '-', 'LineWidth', 1);
    % plot(v_ax, vr_numer, r_numer / R, 'Color', colors(2, :), 'Linestyle', 'None', 'Marker', 'o', 'MarkerSize', 3);
    plot(v_ax, vr_interp, y_interp / R, 'Color', colors(2, :), 'Linestyle', '--', 'LineWidth', 1);
    plot(v_ax, [0, 0], [-1, 1], 'Color', [0.75, 0.75, 0.75] );               % grey
        
    w_ax.XLim  = w_lim;  v_ax.XLim  = v_lim;
    w_ax.YLim = [-1, 1]; v_ax.YLim = [-1, 1];
    grid(w_ax, 'minor'); grid(v_ax, 'minor');
    
    if ii == 1
        legend(w_ax, 'Analytical', 'Numerical', 'Location', 'northeast', 'FontSize', 6);
        legend(v_ax, 'Analytical', 'Numerical', 'Location', 'northwest', 'FontSize', 6);
    end
    
    for jj = 1 : length(ax_all)
        if ii == 1
            ylabel(ax_all{jj}, '$t$ = $0$', 'interpreter', 'latex');
        elseif ii == 2
            ylabel(ax_all{jj}, ['$t$ = $T$ / ', num2str(t_steps)], 'interpreter', 'latex');
        elseif ii == t_steps + 1
            ylabel(ax_all{jj}, '$t$ = $T$', 'interpreter', 'latex');
        else
            ylabel(ax_all{jj}, ['$t$ = ', num2str(ii-1), ' $T$ / ', num2str(t_steps)], 'interpreter', 'latex');
        end
    end
end

sgtitle(w_fig, ['Axial Velocity $w(r, z=',  num2str(z_coord), ', t)$'], 'interpreter', 'latex');
sgtitle(v_fig, ['Radial Velocity $v(r, z=', num2str(z_coord), ', t)$'], 'interpreter', 'latex');

saveas(w_fig, [sim_dir, '/exact-numer_axial-velo-profiles.png'])
saveas(v_fig, [sim_dir, '/exact-numer_radial-velo-profiles.png'])