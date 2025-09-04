clear; clc; close all;
% -----------------------------------------------
% I. Define Parameters
% -----------------------------------------------
% A. Geometric Parameters 
L = 0.450;       % m, Pipe length
D_in = 0.016;    % m, Pipe inner diameter
D_out = 0.019;   % m, Pipe outer diameter

% B. Material Properties 
% SS316
rho_s = 8000;   % kg/m^3, Density of solid
cs = 500;       % J/(kg K), Specific heat of solid
ks = 16.3;      % W/(m K), Thermal conductivity of solid
% Water
rho_f = 997;    % kg/m^3, Density of fluid
cf = 4180;      % J/(kg K), Specific heat of fluid
mu = 0.7972*10^-3; % Pa.s, Viscosity of water
k = 0.6;        % W/(m K), Thermal conductivity of water

% C. Simulation Parameters
Nx = 50;
dx = L / (Nx - 1);
t_final = 30;
dt = 0.1;
Nt = ceil(t_final / dt);
tolerance = 1e-6;
max_iter = 100;

% D. Power Levels, Inlet Temps, and Flow Rates
power_levels = 100:100:800;
num_cases = length(power_levels);

% LPM  corresponding to power levels ---
lpm_values = [4.3, 4.3, 4.2, 4.2, 4.2, 4.2, 4.1, 4.1]; % L.P.M, Flow rate

% Inlet temperatures in Celsius for each power level
inlet_temps_C = [32.77, 32.33, 31.87, 31.49, 31.13, 30.78, 30.41, 30.14];

% Experimental outlet temperatures
Outlet_experimental = [33.02, 32.85, 32.66, 32.54, 32.39, 32.34, 32.22, 32.23];

% E. Derived Parameters (Constant ones)
P_in = pi * D_in;                               % m, Inner perimeter
P_out = pi * D_out;                             % m, Outer perimeter
As = pi/4 * (D_out^2 - D_in^2);                 % m^2, Solid cross-sectional area
Af = pi/4 * D_in^2;                             % m^2, Fluid cross-sectional area

% Insulating material and ambient condition
k_insul = 0.22;         % W/(m K), Insulation thermal conductivity
R2 = D_out / 2;         % m, Insulation inner radius = pipe outer radius
R3 = 116 / 1000;        % m, Insulation outer radius
ha = 5;                 % W/(m^2 K), Assumed ambient convection coeff (natural air)
Ta = 25 + 273.15;       % K, Assumed ambient temperature (25 C)
R_loss = log(R3 / R2) / (2 * pi * k_insul) + 1 / (ha * 2 * pi * R3); % Thermal resistance of insulation + ambient

% F. FDM Coefficients (Constant ones)
Cs1 = rho_s * cs / dt;
Cs2 = ks / dx^2;
Cs4 = 1 / (R_loss * As);


% Initialize storage arrays
x_nodes = linspace(0, L, Nx);                 % Spatial grid nodes
time_vec = linspace(0, t_final, Nt + 1);      % Time vector for plotting
Ts_all = zeros(Nx, Nt + 1, num_cases);
Tf_all = zeros(Nx, Nt + 1, num_cases);
Tf_out_all = zeros(Nt + 1, num_cases);
error_results = zeros(num_cases, 4);

% -----------------------------------------------
% II. Simulation Loop
% -----------------------------------------------
% Loop over different power levels
for p = 1:num_cases
    power = power_levels(p);
    
    % --- Set inlet temperature and LPM for the current power level ---
    T_f_inlet_C = inlet_temps_C(p);
    T_f_inlet = T_f_inlet_C + 273.15;
    T_initial = T_f_inlet; % Assume initial pipe/fluid temp matches the new inlet
    LPM = lpm_values(p); % Get the LPM for the current case
    
    fprintf('Running simulation for Power = %d W, Inlet Temp = %.2f C, Flow Rate = %.1f LPM\n', power, T_f_inlet_C, LPM);
    
    % ---  All flow-dependent parameters are now calculated inside the loop ---
    mass_flow_rate = (LPM / 60000) * rho_f;        % Convert LPM to kg/s
    u = mass_flow_rate / (rho_f * Af);            % Velocity of fluid
    Re = (rho_f * u * D_in) / mu;                 % Reynolds number
    Pr = (mu * cf) / k;                           % Prandtl number of fluid
    f = (0.790 * log(Re) - 1.64)^(-2);             % Friction factor for smooth pipe
    L_td = 10 * D_in;                             % Thermal entrance length
    Nu_td = 0.036 * (Re^0.8) * (Pr^0.385) * (L / D_in)^(0.054); % Ghajar and Tam relation
    Nu_fd = ((f/8) * (Re - 1000) * Pr) / (1 + 12.7 * sqrt(f/8) * (Pr^(2/3) - 1)); % Gnielinski correlation
    Nu = (Nu_td * L_td + Nu_fd * (L - L_td)) / L; % Average Nusselt number
    hsf = (Nu * k) / D_in;                        % Convective heat transfer coefficient
    %disp(hsf)
    %disp(Nu)


    % FDM Coefficients that depend on flow
    Cf1 = rho_f * cf / dt;
    Cf2 = rho_f * cf * u / dx;
    Cf3 = hsf * P_in / Af;
    Cs3 = hsf * P_in / As;

    % Volumetric heat source based on input power
    q_prime_prime_in = power / (pi * D_out * L);
    Q_source_s_val = (q_prime_prime_in * P_out) / As;
    
    % Initialize temperature fields for this case at t=0;
    Ts = ones(Nx, 1) * T_initial;      % Pipe wall temperature array (K)
    Tf = ones(Nx, 1) * T_initial;      % Fluid temperature array (K)
    Tf(1) = T_f_inlet;                 % Set fluid inlet temperature (K)
    
    % Store initial conditions
    Ts_all(:, 1, p) = Ts;
    Tf_all(:, 1, p) = Tf;
    Tf_outlet_vs_time = zeros(Nt + 1, 1);
    Tf_outlet_vs_time(1) = Tf(end);
    
    % Time marching loop with Picard iteration
    for n = 1:Nt
        Ts_old = Ts_all(:, n, p);
        Tf_old = Tf_all(:, n, p);
        
        % Iterative solver for coupled equations at time t+dt
        Ts_iter = Ts_old;
        Tf_iter = Tf_old;
        
        for iter = 1:max_iter
            Ts_prev = Ts_iter;
            Tf_prev = Tf_iter;
            
            % --- Solve for Ts_iter ---
            diag_s = zeros(Nx, 1);
            sub_diag_s = zeros(Nx - 1, 1);
            sup_diag_s = zeros(Nx - 1, 1);
            Bs_rhs = zeros(Nx, 1);

            % Boundary node i=1 (x=0), insulated
            diag_s(1) = Cs1 + 2 * Cs2 + Cs3 + Cs4;
            sup_diag_s(1) = -2 * Cs2;
            Bs_rhs(1) = Cs1 * Ts_old(1) + Q_source_s_val + Cs3 * Tf_prev(1) + Ta * Cs4;
            
            % Internal nodes i=2 to Nx-1
            for i = 2:Nx-1
                diag_s(i) = Cs1 + 2 * Cs2 + Cs3 + Cs4;
                sub_diag_s(i-1) = -Cs2;
                sup_diag_s(i) = -Cs2;
                Bs_rhs(i) = Cs1 * Ts_old(i) + Q_source_s_val + Cs3 * Tf_prev(i) + Ta * Cs4;
            end
            
            % Boundary node i=Nx (x=L), insulated
            diag_s(Nx) = Cs1 + 2 * Cs2 + Cs3 + Cs4;
            sub_diag_s(Nx-1) = -2 * Cs2;
            Bs_rhs(Nx) = Cs1 * Ts_old(Nx) + Q_source_s_val + Cs3 * Tf_prev(Nx) + Ta * Cs4;            
            As_matrix = diag(diag_s) + diag(sub_diag_s, -1) + diag(sup_diag_s, 1);
            Ts_iter = As_matrix \ Bs_rhs;
            
            % --- Solve for Tf_iter ---
            Tf_iter(1) = T_f_inlet; % Dirichlet BC
            num_unknowns_f = Nx - 1;

            diag_f = zeros(num_unknowns_f, 1);
            sub_diag_f = zeros(num_unknowns_f - 1, 1);
            Bf_rhs = zeros(num_unknowns_f, 1);

            % Node i=2 (first unknown)
            diag_f(1) = Cf1 + Cf2 + Cf3;
            Bf_rhs(1) = Cf1 * Tf_old(2) + Cf3 * Ts_iter(2) + Cf2 * T_f_inlet;


            % Internal nodes i=3 to Nx
            for i = 2:num_unknowns_f
                diag_f(i) = Cf1 + Cf2 + Cf3;
                sub_diag_f(i-1) = -Cf2;
                Bf_rhs(i) = Cf1 * Tf_old(i+1) + Cf3 * Ts_iter(i+1);
            end
            Af_matrix_internal = diag(diag_f) + diag(sub_diag_f, -1);
            Tf_iter(2:Nx) = Af_matrix_internal \ Bf_rhs;
            
            % Convergence check
            if norm(Ts_iter - Ts_prev, inf) < tolerance && norm(Tf_iter - Tf_prev, inf) < tolerance
                break;
            end
        end
        
        % Store results for current time step
        Ts_all(:, n + 1, p) = Ts_iter;
        Tf_all(:, n + 1, p) = Tf_iter;
        Tf_outlet_vs_time(n + 1) = Tf_iter(end);
    end
    
    % Store outlet temperature vs. time for this power level (in Celsius)
    Tf_out_all(:, p) = Tf_outlet_vs_time - 273.15;
end

% -----------------------------------------------
% III. Plotting and Reporting Results
% -----------------------------------------------
% 1. Summary Plot: Water Outlet Temperature vs. Time for All Power Levels
figure;
hold on;
colors = lines(num_cases); % Generate distinct colors for plotting
for p = 1:num_cases
    plot(time_vec, Tf_out_all(:, p), 'LineWidth', 1.5, ...
        'Color', colors(p,:), ...
        'DisplayName', sprintf('%d W (Inlet %.2f C)', power_levels(p), inlet_temps_C(p)));
end
hold off;
xlabel('Time (s)');
ylabel('Outlet Water Temperature (°C)');
title('Water Outlet Temperature vs. Time for All Power Levels');
legend('show', 'Location', 'southeast');
grid on;

% 2. Compare Simulated vs. Experimental Outlet Temperature for each power level
for p = 1:num_cases
    power = power_levels(p);
    exp_temp = Outlet_experimental(p);
    sim_temp_final = Tf_out_all(end, p);
    percentage_error = abs((sim_temp_final - exp_temp) / exp_temp) * 100;
    error_results(p, :) = [power, exp_temp, sim_temp_final, percentage_error]; 

    figure;
    hold on;
    plot(time_vec, Tf_out_all(:, p), 'b-', 'LineWidth', 2, 'DisplayName', 'Simulated Outlet Temp');
    yline(exp_temp, 'r--', 'LineWidth', 2, 'DisplayName', sprintf('Experimental Steady State (%.2f °C)', exp_temp));
    hold off;
    xlabel('Time (s)');
    ylabel('Water  Temperature (°C)');
    title(sprintf('Water temp vs time at %d W | Error: %.2f%%', power, percentage_error));
    legend('show', 'Location', 'southeast');
    grid on;
end

% 3. Display Error Summary Table
fprintf('\n--- Comparison of Simulated vs. Experimental Outlet Temperatures ---\n');
fprintf('Power (W) | Exp. Temp (°C) | Sim. Temp (°C) | Error (%%)\n');
fprintf('------------------------------------------------------------------\n');
for i = 1:num_cases
    fprintf('%-9d | %-14.2f | %-14.2f | %.2f\n', ...
        error_results(i, 1), error_results(i, 2), error_results(i, 3), error_results(i, 4));
end
fprintf('------------------------------------------------------------------\n');
