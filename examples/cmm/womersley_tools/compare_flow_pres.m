function compare_flow_pres(sim_dir, num_cyc, z_in, z_out, inlet_data, outlet_data, p0, mu, rho, R, c_n, B_n, Q_n, G_n, g_n, T, n_modes, sol_idx)

conversion = 1333.2;

colors = [     0, 0.4470, 0.7410; 0.8500, 0.3250, 0.0980; ...
          0.4940, 0.1840, 0.5560; 0.4660, 0.6740, 0.1880; ...
          0.3010, 0.7450, 0.9330; 0.6350, 0.0780, 0.1840];
      
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

t_numer = inlet_data(sol_idx, 2);

q_in_numer  = -inlet_data(sol_idx, 4); % flip sign for unit normal
q_out_numer = outlet_data(sol_idx, 4);

p_in_numer  =  inlet_data(sol_idx, 5);
p_out_numer = outlet_data(sol_idx, 5);


% Plot analytical vs. numerical flows
figure; hold on;
plot(t_exact, real(q_in_exact),  'Color', colors(1, :), 'Linestyle', '-');
plot(t_exact, real(q_out_exact), 'Color', colors(2, :), 'Linestyle', '-');
plot(t_numer,  q_in_numer, 'Color', colors(1, :), 'Linestyle', '--');
plot(t_numer, q_out_numer, 'Color', colors(2, :), 'Linestyle', '--');

xlabel('Time (s)'); ylabel('Flow (ml/s)'); xlim([0, T]);
legend('Inlet / Analytical', 'Outlet / Analytical', 'Inlet / Numerical', 'Outlet / Numerical');
saveas(gcf, [sim_dir, '/exact-numer_cap-flows.png'])

% Plot analytical vs. numerical pressures
figure; hold on;
plot(t_exact, real(p_in_exact)  / conversion, 'Color', colors(1, :), 'Linestyle', '-');
plot(t_exact, real(p_out_exact) / conversion, 'Color', colors(2, :), 'Linestyle', '-');
plot(t_numer,  p_in_numer / conversion, 'Color', colors(1, :), 'Linestyle', '--');
plot(t_numer, p_out_numer / conversion, 'Color', colors(2, :), 'Linestyle', '--');

xlabel('Time (s)'); ylabel('Pressure (mm Hg)'); xlim([0, T]);
legend('Inlet / Analytical', 'Outlet / Analytical', 'Inlet / Numerical', 'Outlet / Numerical');
saveas(gcf, [sim_dir, '/exact-numer_cap-pressures.png'])



