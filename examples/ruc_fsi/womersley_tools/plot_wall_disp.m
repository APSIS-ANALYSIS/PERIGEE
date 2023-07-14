function plot_wall_disp(rho, R, c_n, g_n, B_n, G_n, T, n_modes, t_steps)

colors = [     0, 0.4470, 0.7410; 0.8500, 0.3250, 0.0980; ...
          0.4940, 0.1840, 0.5560; 0.4660, 0.6740, 0.1880; ...
          0.3010, 0.7450, 0.9330; 0.6350, 0.0780, 0.1840];

dt = T / t_steps;

omega = 2 * pi / T;                                 % base angular frequency

% Only entries starting at index 2 are meaningful
c_R = 1 ./ real(1 ./ c_n);                          % real wave speed (dispersion coefficient)
L_n = c_R * 2 * pi ./ ((0 : n_modes - 1) * omega);  % wavelengths

z = linspace(0, L_n(2), 100);                       % axial coordinate

xi_fig  = figure;  xi_ax = axes;  hold on;
eta_fig = figure; eta_ax = axes;  hold on;

xi_lim  = [-0.8, 0.8]; eta_lim = [-1.6e-3, 1.6e-3];

t_labs = cell(1, t_steps + 1);
t_labs{1} = '$t$ = $0$';
t_labs{2} = ['$t$ = $T$ / ', num2str(t_steps)];
for ii = 3 : t_steps
    t_labs{ii} = ['$t$ = ', num2str(ii-1), ' $T$ / ', num2str(t_steps)];
end
t_labs{t_steps + 1} = '$t$ = $T$';

for ii = 1 : (t_steps + 1)
    
    t = (ii - 1) * dt;
    
    xi  = zeros(1, length(z));                      % axial displacement
    eta = zeros(1, length(z));                      % radial displacement
    
    for k = 2 : n_modes
        
        n = k - 1;
        
        xi  = xi + 1j * B_n(k) / ( rho * c_n(k) * n * omega ) * ( G_n(k) - 1 ) * ...
              exp(1j * n * omega * (t - z / c_n(k)) );
        eta = eta + B_n(k) * R / ( 2 * rho * c_n(k)^2 ) * ( 1 - G_n(k) * g_n(k) ) * ...
              exp(1j * n * omega * (t - z / c_n(k)) ); 
    end
    
    plot(xi_ax,  z / L_n(2), real(xi),  'Color', colors(ii, :), 'Linestyle', '-');
    plot(eta_ax, z / L_n(2), real(eta), 'Color', colors(ii, :), 'Linestyle', '-');
    
end

xi_ax.YLim = xi_lim; eta_ax.YLim = eta_lim;

xlabel(xi_ax,  'z / $\lambda$', 'interpreter', 'latex');
xlabel(eta_ax, 'z / $\lambda$', 'interpreter', 'latex');
title(xi_ax,  'Axial Wall Disp $\xi(z, t)$',   'interpreter', 'latex');
title(eta_ax, 'Radial Wall Disp $\eta(z, t)$', 'interpreter', 'latex');

legend(xi_ax,  t_labs, 'interpreter', 'latex', 'NumColumns', 3);
legend(eta_ax, t_labs, 'interpreter', 'latex', 'NumColumns', 3);

saveas(xi_fig,  'exact_axial-wall-disp.png');
saveas(eta_fig, 'exact_radial-wall-disp.png');