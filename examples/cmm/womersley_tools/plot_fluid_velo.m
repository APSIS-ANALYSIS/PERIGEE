function plot_fluid_velo(mu, rho, R, c_n, B_n, Q_n, G_n, T, n_modes, z_steps, t_steps)

colors = [     0, 0.4470, 0.7410; 0.8500, 0.3250, 0.0980; ...
          0.4940, 0.1840, 0.5560; 0.4660, 0.6740, 0.1880; ...
          0.3010, 0.7450, 0.9330; 0.6350, 0.0780, 0.1840];

dt = T / t_steps;

omega = 2 * pi / T;                                 % base angular frequency

% Only entries starting at index 2 are meaningful
c_R = 1 ./ real(1 ./ c_n);                          % real wave speed (dispersion coefficient)
L_n = c_R * 2 * pi ./ ((0 : n_modes - 1) * omega);  % wavelengths
assignin('base', 'L_n', L_n);

dz = L_n(2) / z_steps;

x = -R : 0.01 : R;  % radius vector
r = abs(x);

w_fig = figure;    v_fig = figure;
w_lim = [-10, 180]; v_lim = [-0.01, 0.11];

% Define offset down the tube length
w_offset = (0 : z_steps) * 30;
v_offset = (0 : z_steps) * 0.02;

% Define xticks
w_xticks = sort( [w_offset + 10,    w_offset + 20  ] );
v_xticks = sort( [v_offset + 0.01] );

% Define xticklabels
w_xticklabs = {'10', '20'};
v_xticklabs = {'0.01'};

% Define labels for z steps
w_zsteplabs = repmat({''}, 1, ( length(w_xticklabs) + 1 ) * (z_steps + 1) );
v_zsteplabs = repmat({''}, 1, ( length(v_xticklabs) + 1 ) * (z_steps + 1) );

for ii = 1 : (z_steps + 1)
    w_idx = (length(w_xticklabs) + 1) * (ii - 1) + 1;
    v_idx = (length(v_xticklabs) + 1) * (ii - 1) + 1;
    
    if ii == 1
        w_zsteplabs{w_idx} = '$z$ = $0$';
        v_zsteplabs{v_idx} = '$z$ = $0$';
    elseif ii == 2
        w_zsteplabs{w_idx} = ['$z$ = $\lambda$ / ', num2str(z_steps)];
        v_zsteplabs{v_idx} = ['$z$ = $\lambda$ / ', num2str(z_steps)];
    elseif ii == z_steps + 1
        w_zsteplabs{w_idx} = '$z$ = $\lambda$';
        v_zsteplabs{v_idx} = '$z$ = $\lambda$';
    else
        w_zsteplabs{w_idx} = ['$z$ = ', num2str(ii-1), ' $\lambda$ / ', num2str(z_steps)];
        v_zsteplabs{v_idx} = ['$z$ = ', num2str(ii-1), ' $\lambda$ / ', num2str(z_steps)];
    end
end

w_xticklabs = repmat( w_xticklabs, 1, z_steps + 1);
v_xticklabs = repmat( v_xticklabs, 1, z_steps + 1);

for ii = 1 : (t_steps + 1)
    
    t = (ii - 1) * dt;
    w_ax = subplot(t_steps + 1, 1, ii, 'Parent', w_fig); hold(w_ax, 'on');
    v_ax = subplot(t_steps + 1, 1, ii, 'Parent', v_fig); hold(v_ax, 'on');
    
    ax_all = {w_ax, v_ax};
    
    for jj = 1 : (z_steps + 1)
        
        z = dz * (jj - 1);
        
        % Initialize with Poiseuille solution (0th mode)
        w_rigid = 2 * Q_n(1) * (R^2 - r.^2) / (pi * R^4);   % axial velocity (rigid theory)
        w = 2 * Q_n(1) * (R^2 - r.^2) / (pi * R^4);         % axial velocity
        
        v = zeros(1, length(r));                            % radial velocity
        
        for k = 2 : n_modes
            n = k - 1;
            Omega_n  = R * sqrt(rho * n * omega / mu);   % Womersley number
            Lambda_n = 1j^1.5 * Omega_n;
            
            w_rigid = w_rigid + Q_n(k) * Lambda_n  / (pi * R^2) * ...
                  (besselj(0, Lambda_n) - besselj(0, Lambda_n * r / R)) / ...
                  (besselj(0, Lambda_n) * Lambda_n - 2 * besselj(1, Lambda_n)) * ...
                  exp(1j * n * omega * t);
          
            w = w + B_n(k) / ( rho * c_n(k) ) * ...
                ( 1 - G_n(k) * besselj(0, Lambda_n * r / R) / besselj(0, Lambda_n) ) * ...
                exp(1j * n * omega * (t - z / c_n(k)) );
            v = v + B_n(k) * 1j * n * omega * R / ( 2 * rho * c_n(k)^2 ) .* ...
                ( r / R - G_n(k) * 2 * besselj(1, Lambda_n * r / R) / (Lambda_n * besselj(0, Lambda_n)) ) * ...
                exp(1j * n * omega * (t - z / c_n(k)) );
        end

        plot(w_ax, w_offset(jj) + real(w_rigid), x / R,   'Color', colors(jj, :), 'Linestyle', '--');
        plot(w_ax, w_offset(jj) + real(w),       x / R,   'Color', colors(jj, :), 'Linestyle', '-' );
        plot(w_ax, [w_offset(jj), w_offset(jj)], [-1, 1], 'Color', [0.75, 0.75, 0.75] ); % grey
        
        plot(v_ax, v_offset(jj) + real(v),       x / R, 'Color', colors(jj, :), 'Linestyle', '-' );
        plot(v_ax, [v_offset(jj), v_offset(jj)], [-1, 1], 'Color', [0.75, 0.75, 0.75] ); % grey

    end
    
    w_ax.XLim  = w_lim;    v_ax.XLim  = v_lim;
    w_ax.XTick = w_xticks; v_ax.XTick = v_xticks;
    
    if ii == 1
        w_ax.XTickLabel = w_xticklabs; v_ax.XTickLabel = v_xticklabs;
    elseif ii == t_steps + 1
        w_ax.XTick = sort( [w_xticks, w_offset] );
        v_ax.XTick = sort( [v_xticks, v_offset] );
        
        w_ax.TickLabelInterpreter = 'latex'; v_ax.TickLabelInterpreter = 'latex';
        w_ax.XTickLabel = w_zsteplabs;       v_ax.XTickLabel = v_zsteplabs;
    else
        w_ax.XTickLabel = []; v_ax.XTickLabel = [];
    end
    
    grid(w_ax, 'minor'); grid(v_ax, 'minor');
    
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

sgtitle(w_fig, 'Axial Velocity $w(r, z, t)$',  'interpreter', 'latex');
sgtitle(v_fig, 'Radial Velocity $v(r, z, t)$', 'interpreter', 'latex');
