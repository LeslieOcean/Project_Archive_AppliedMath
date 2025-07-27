clear all; clc;

itera = 15;
domega = [0; 0]; %w1, w2, w3
k = 1; %linspace(0, 2, itera);
%IC = pi*rand(itera,2);
tspan = linspace(0,100,5000);
phi_12_list = zeros(itera,1);
phi_23_list = zeros(itera,1);

% for i=1:itera %for different k
%     phi_12 = 0;
%     phi_23 = 0;
% for j=1:itera %for different ICs
phi0 = [pi/2; 0];
[t,phi] = ode45(@(t,phi)kuramoto_3oscillator_1(t,phi,domega,k), tspan, phi0);
% phi_12 = phi_12 + phi(end,1);
% phi_23 = phi_23 + phi(end,2);
% end
% phi_12_list(i) = phi_23/itera;
% phi_23_list(i) = phi_23/itera;
% end

figure(1);
plot(t,phi(:,1),color='blue');
hold on;
plot(t,phi(:,2),color='red');
% legend('sin(\theta_1)','sin(\theta_2)')

% figure(2);
% plot(t,cos(phi));
% legend('cos(\phi)')
% 
% figure(3);
% plot(t,phi);
% legend('\phi');

%%
% Parameters (example values - replace with your actual parameters)
clear all; clc;
domega = [0; 0];   % Vector [domega(1); domega(2)]
k = 1;                  % Coupling strength

% Define the system of equations
fun1 = @(phi) [
    (-2*sin(phi(1)) - sin(phi(1) + phi(2)) + sin(phi(2)))*k/3 + domega(1);
    (sin(phi(1)) - sin(phi(1) + phi(2)) - 2*sin(phi(2)))*k/3 + domega(2)
];

fun2 = @(phi) [
    (-5*sin(phi(1)) + 2*sin(phi(2)))*k/6 + domega(1);
    (2*sin(phi(1)) - 5*sin(phi(2)))*k/4 + domega(2)
];

% Settings for fsolve
options = optimoptions('fsolve', 'Display', 'off');  % Suppress output

% Grid of initial guesses in [0, 2*pi]
[x0, y0] = meshgrid(linspace(0, 2*pi, 40), linspace(0, 2*pi, 40));
solutions = [];

% Solve for each initial guess
for i = 1:numel(x0)
    phi0 = [x0(i); y0(i)];                % Initial guess
    [phi_sol, ~, exitflag] = fsolve(fun2, phi0, options);
    
    % Check if solution is valid and unique
    if exitflag > 0
        phi_sol = mod(phi_sol, 2*pi);     % Wrap to [0, 2*pi)
        % Check for uniqueness (within tolerance)
        if isempty(solutions) || all(min(vecnorm(solutions - phi_sol, 2, 1)) > 1e-4)
            solutions = [solutions, phi_sol];
        end
    end
end

% Display all unique fixed points
disp('Fixed points found (phi1, phi2):');
disp(solutions');


%%
clear all; clc;

k_range = linspace(0, 1, 50);         % 50 values of k
w_range = linspace(-0.5, 0.5, 50);       % 50 values of delta_omega
n_trials = 40;                          % Number of initial conditions
tspan = [0 100];                          % Time span for ODE
init_cond = 2 * pi * rand(3, n_trials);  % Initial conditions

% === Output matrix ===
phi12_avg = zeros(length(k_range), length(w_range));  % for contour plot
phi23_avg = zeros(length(k_range), length(w_range));
dphi_avg = zeros(length(k_range), length(w_range));

% === Nested loops over k and delta_omega ===
for ki = 1:length(k_range)
    k = k_range(ki);

    for wi = 1:length(w_range)
        delta_omega = w_range(wi);
        omega1 = 1;
        omega2 = omega1 - delta_omega;
        omega3 = omega2 - delta_omega;
        phi12dot_vals = zeros(n_trials, 1);
        phi23dot_vals = zeros(n_trials, 1);
        dphi_dot_vals = zeros(n_trials, 1);
        for j = 1:n_trials
            theta0_rand = init_cond(:, j);

            % ODE system
            ode_fun2 = @(t, theta) [
                omega1 + k/2 * (sin((theta(2)-theta(1))));
                omega2 + k/3 * (sin((theta(3)-theta(2))) + sin((theta(1)-theta(2))));
                omega3 + k/2 * (sin((theta(2)-theta(3))));
            ];

            ode_fun1 = @(t, theta) [
                omega1 + k/3 * (sin(theta(2)-theta(1))+sin(theta(3)-theta(1)));
                omega2 + k/3 * (sin(theta(1)-theta(2))+sin((theta(3)-theta(2))));
                omega3 + k/3 * (sin(theta(2)-theta(3))+sin(theta(1)-theta(3)));
            ];


            % Solve ODE
            [~, theta] = ode45(ode_fun2, tspan, theta0_rand);
            phi12 = theta(end,1) - theta(end,2);
            phi23 = theta(end,2) - theta(end,3);

            % Phase difference at final time
            dphi12 = delta_omega + k/6 * (-5*sin(phi12) + 2*sin(phi23)); %config2
            dphi23 = delta_omega + k/6 * (2*sin(phi12) - 5*sin(phi23));
            %dphi12 = delta_omega + k/3 * (-2*sin(phi12)+sin(phi23)-sin(phi12+phi23))
            %dphi23 = delta_omega + k/3 * (-2*sin(phi23)+sin(phi12)-sin(phi12+phi23))
            ddphi = dphi12 + dphi23;
            phi12dot_vals(j) = dphi12;
            phi23dot_vals(j) = dphi23;
            dphi_dot_vals(j) = ddphi;
        end

        % Average over all trials
        phi12_avg(ki, wi) = mean(phi12dot_vals);
        phi23_avg(ki, wi) = mean(phi23dot_vals);
        dphi_avg(ki, wi) = mean(phi23dot_vals);
    end
end

% === Contour plot ===
figure;
contourf(w_range, k_range, phi12_avg, 20, 'LineColor', 'none');
colorbar;
% colormap(parula);      % Better than default 'jet'
colormap(turbo);       % High contrast and smooth — recommended
% colormap('hot');       % Strong warm contrast
% colormap('viridis');   % If available
xlabel('\Delta\omega');
ylabel('Coupling Strength k [rad/s]');

% title('Contour Plot of $\dot{\phi}$ vs k and \Delta\omega', 'Interpreter', 'latex');
% === Contour plot ===
figure;
contourf(w_range, k_range, phi23_avg, 20, 'LineColor', 'none');
colorbar;
% colormap(parula);      % Better than default 'jet'
colormap(turbo);       % High contrast and smooth — recommended
% colormap('hot');       % Strong warm contrast
% colormap('viridis');   % If available
xlabel('\Delta\omega');
ylabel('Coupling Strength k [rad/s]');
% title('Contour Plot of $\dot{\phi}$ vs k and \Delta\omega', 'Interpreter', 'latex');


% title('Contour Plot of $\dot{\phi}$ vs k and \Delta\omega', 'Interpreter', 'latex');
% === Contour plot ===
figure;
contourf(w_range, k_range, dphi_avg, 20, 'LineColor', 'none');
colorbar;
% colormap(parula);      % Better than default 'jet'
colormap(turbo);       % High contrast and smooth — recommended
% colormap('hot');       % Strong warm contrast
% colormap('viridis');   % If available
xlabel('\Delta\omega');
ylabel('Coupling Strength k [rad/s]');
% title('Contour Plot of $\dot{\phi}$ vs k and \Delta\omega', 'Interpreter', 'latex');