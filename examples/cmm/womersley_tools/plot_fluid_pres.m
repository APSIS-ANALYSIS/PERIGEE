function plot_fluid_pres(p0, mu, R, c_n, B_n, Q_n, T, n_modes, t_steps)

conversion = 1333.2;

colors = [     0, 0.4470, 0.7410; 0.8500, 0.3250, 0.0980; ...
          0.4940, 0.1840, 0.5560; 0.4660, 0.6740, 0.1880; ...
          0.3010, 0.7450, 0.9330; 0.6350, 0.0780, 0.1840];

dt = T / t_steps;

omega = 2 * pi / T;                                 % base angular frequency

% Only entries starting at index 2 are meaningful
c_R = 1 ./ real(1 ./ c_n);                          % real wave speed (dispersion coefficient)
L_n = c_R * 2 * pi ./ ((0 : n_modes - 1) * omega);  % wavelengths

z = linspace(0, L_n(2), 100);                       % axial coordinate

figure; hold on;

t_labs = cell(1, t_steps + 1);
t_labs{1} = '$t$ = $0$';
t_labs{2} = ['$t$ = $T$ / ', num2str(t_steps)];
for ii = 3 : t_steps
    t_labs{ii} = ['$t$ = ', num2str(ii-1), ' $T$ / ', num2str(t_steps)];
end
t_labs{t_steps + 1} = '$t$ = $T$';

for ii = 1 : (t_steps + 1)
    
    t = (ii - 1) * dt;
    
    % Initialize with Poiseuille solution (0th mode)
    p = p0 -8 * mu * Q_n(1) / (pi * R^4) * z;
    
    for k = 2 : n_modes
        
        n = k - 1;
        
        p = p + B_n(k) * exp(1j * n * omega * (t - z / c_n(k)) );
    end
    
    plot(z / L_n(2), real(p) / conversion, 'Color', colors(ii, :), 'Linestyle', '-');
    
end

xlabel('z / $\lambda$', 'interpreter', 'latex');
ylabel('mm Hg', 'interpreter', 'latex')
title('Fluid Pressure $P(r, z, t)$', 'interpreter', 'latex');

legend(t_labs, 'interpreter', 'latex');

