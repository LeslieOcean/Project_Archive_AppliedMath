clear all; clc;

syms h2 u v;

rho=1.76;
k=5;
theta=1.33;
delta=1.5;

a=0.02; %(0,1)
h1=0.24;
theta1=0.21; %(0,1)

A = - (rho * (theta1 - a)) / (theta * k);

B = (rho * (theta1 - a)) / theta * (1 - delta / k) + (a * delta * rho) / (theta * k);

C = (theta1 - a) * (h1 + (rho * delta) / theta) - (a * delta * rho) / theta * (1 - delta / k) - h1;

D = - a * delta * h1 - (a * delta^2 * rho) / theta - h1 * delta;

Delta = 18*A*B*C*D - 4*B^3*D + B^2*C^2 - 4*A*C^3 - 27*A^2*D^2;
fprintf('Discriminant = %.6f\n', Delta);

cubic_u = A*u^3 + B*u^2 + C*u + D;
sol_u = solve(cubic_u == 0, u);
u_cubic_vals = vpa(sol_u, 6);

C1 = (theta1 - a) * (h2 + (rho * delta) / theta) - (a * delta * rho) / theta * (1 - delta / k) - h2;
D1 = - a * delta * h2 - (a * delta^2 * rho) / theta - h2 * delta;
eqns = 18*A*B*C1*D1 - 4*B^3*D1 + B^2*C1^2 - 4*A*C1^3 - 27*A^2*D1^2==0;
S = solve(eqns, h2, MaxDegree=3);
S_vals = vpa(S, 6);
h2_val = double(S_vals(real(S_vals) == S_vals & S_vals > 0));
h2 = h2_val

asym = a*delta/(theta1-a);
point1x = k;
point1y = 0;
point2x = 0;
point2y = rho*delta/theta;

u1 = linspace(0,k,100);
u2 = linspace(asym+0.1, k+1, 100);
v_prey=rho.*(k-u1).*(delta+u1)./(k*theta);
v_predator1=h1.*(delta+u2)./(-a.*(delta+u2)+theta1.*u2)-h1;
v_predator2=h2.*(delta+u2)./(-a.*(delta+u2)+theta1.*u2)-h2;

fp1x = double(u_cubic_vals(u_cubic_vals>0));
fp1y = rho.*(k-fp1x).*(delta+fp1x)./(k*theta);

%% Figure 1
figure(1)
plot(u1, v_prey, 'LineWidth', 1.5);
hold on;
plot(u2, v_predator1, 'LineWidth', 1.5);
plot(point1x, point1y, 'm*', 'MarkerSize', 5);
plot(point2x, point2y, 'b*', 'MarkerSize', 5);
plot(fp1x, fp1y, 'o', 'MarkerSize', 5);

ylim([-1, 5]);
xlim([-1, k + 1]);
xticks([]);
yticks([]);
xlabel('u', 'FontSize', 16);
ylabel('v', 'FontSize', 16);
yline(0, 'k');
xline(0, 'k');
xline(asym, '--b', 'a\delta/(\theta_1-a)', 'FontSize', 14, 'LabelVerticalAlignment', 'bottom', 'LineWidth', 1);
legend({'prey isocline','predator isocline','(k,0)','(0,\rho\delta/\theta)','Equilibria'},'FontSize', 13);
set(gca, 'FontSize', 14);

%% Figure 2
figure(2);
plot(u1, v_prey, 'LineWidth', 1.5);
hold on;
plot(u2, v_predator2, 'LineWidth', 1.5);
plot(point1x, point1y, 'm*', 'MarkerSize', 5);
plot(point2x, point2y, 'b*', 'MarkerSize', 5);
plot(2.63594, 2.58777, 'o', 'MarkerSize', 5);

ylim([-1, 5]);
xlim([-1, k + 1]);
xticks([]);
yticks([]);
xlabel('u', 'FontSize', 16);
ylabel('v', 'FontSize', 16);
yline(0, 'k');
xline(0, 'k');
xline(asym, '--b', 'a\delta/(\theta_1-a)', 'FontSize', 14, 'LabelVerticalAlignment', 'bottom', 'LineWidth', 1);
legend({'prey isocline','predator isocline','(k,0)','(0,\rho\delta/\theta)','Equilibrium'}, 'FontSize', 13);
set(gca, 'FontSize', 14);

%% Figure 3
h3 = 0.5;
v_predator3 = h3 * (delta + u2) ./ (-a .* (delta + u2) + theta1 .* u2) - h3;

figure(3);
plot(u1, v_prey, 'LineWidth', 1.5);
hold on;
plot(u2, v_predator3, 'LineWidth', 1.5);
plot(point1x, point1y, 'm*', 'MarkerSize', 5);
plot(point2x, point2y, 'b*', 'MarkerSize', 5);

ylim([-1, 5]);
xlim([-1, k + 1]);
xticks([]);
yticks([]);
xlabel('u', 'FontSize', 16);
ylabel('v', 'FontSize', 16);yline(0, 'k');
xline(0, 'k');
xline(asym, '--b', 'a\delta/(\theta_1-a)', 'FontSize', 14, 'LabelVerticalAlignment', 'bottom', 'LineWidth', 1);
legend({'prey isocline','predator isocline','(k,0)','(0,\rho\delta/\theta)'}, 'FontSize', 13);
set(gca, 'FontSize', 14);

%% jacobian

syms u v rho k theta1 delta a h theta real 

f = rho * u * (1 - u/k) - (theta * u * v)/(delta + u);
g = (theta1 * u * v)/(delta + u) - a * v - (h * v)/(h + v);

J = jacobian([f; g], [u, v]);
disp('Jacobian matrix:');
pretty(J)

u_star = 2.63594;
v_star = 2.58777;

rho_val = 1.76;
k_val = 5;
theta_val = 1.33;
theta1_val = 0.21;
delta_val = 1.5;
a_val = 0.02;
h_val = 0.24;

J_subs = subs(J, [u, v, rho, k, theta1, delta, a, h, theta], ...
                 [u_star, v_star, rho_val, k_val, theta1_val, delta_val, a_val, h_val, theta_val]);

J_num = double(J_subs);

disp('Numeric Jacobian matrix at equilibrium:');
disp(J_num);

eig_vals = eig(J_num);
disp('Eigenvalues of the Jacobian:');
disp(eig_vals);


%%
syms u v rho k theta1 delta a h theta

f = rho * u * (1 - u/k) - (theta1 * u * v)/(delta + u);
g = (theta1 * u * v)/(delta + u) - a * v - (h * v)/(h + v);
J = jacobian([f; g], [u, v]);

rho_val = 1.76;
k_val = 5;
theta_val = 1.33;
theta1_val = 0.21;
delta_val = 1.5;
a_val = 0.02;

A = - (rho_val * (theta1_val - a_val)) / (theta_val * k_val);
B = (rho_val * (theta1_val - a_val)) / theta_val * (1 - delta_val / k_val) + ...
    (a_val * delta_val * rho_val) / (theta_val * k_val);

h_values = linspace(0.1, 0.33, 10);  % adjust upper limit if needed
arg = zeros(length(h_values), 1);
eigenvalue_storage = zeros(length(h_values), 2); % store dominant eigenvalues

for i = 1:length(h_values)
    h_i = h_values(i);

    C = (theta1_val - a_val) * (h_i + (rho_val * delta_val) / theta_val) - ...
        (a_val * delta_val * rho_val) / theta_val * (1 - delta_val / k_val) - h_i;
    D = - a_val * delta_val * h_i - (a_val * delta_val^2 * rho_val) / theta_val - h_i * delta_val;

    cubic_roots = roots([A, B, C, D]);
    real_pos_roots = real(cubic_roots(abs(imag(cubic_roots)) < 1e-6 & real(cubic_roots) > 0));

    if isempty(real_pos_roots)
        arg(i) = NaN;
        continue
    end

    u_star = real_pos_roots(1);
    v_star = rho_val * (k_val - u_star) * (delta_val + u_star) / (k_val * theta_val);

    J_subs = subs(J, [u, v, rho, k, theta1, delta, a, h, theta], ...
                     [u_star, v_star, rho_val, k_val, theta1_val, delta_val, a_val, h_i, theta_val]);

    J_num = double(J_subs);
    eig_vals = eig(J_num);
    eigenvalue_storage(i, :) = eig_vals.';

    ev = eig_vals(1); % or use the one with max real part
    arg(i) = atan2(imag(ev), real(ev)); % gives phase angle in radians
end

figure;
plot(h_values, arg, 'LineWidth', 1.5);
xlabel('h');
ylabel('arg(\lambda)');
title('Phase of Dominant Eigenvalue vs h');
grid on;

%%
rho = 1.76;
k = 5;
theta = 1.33;
theta1 = 0.21;
delta = 1.5;
a = 0.02;
h = 0.24;

[u_grid, v_grid] = meshgrid(linspace(0, 6, 40), linspace(0, 6, 40));

du = rho * u_grid .* (1 - u_grid / k) - (theta * u_grid .* v_grid) ./ (delta + u_grid);
dv = (theta1 * u_grid .* v_grid) ./ (delta + u_grid) - a * v_grid - (h * v_grid) ./ (h + v_grid);

magnitude = sqrt(du.^2 + dv.^2);
du_norm = du ./ magnitude;
dv_norm = dv ./ magnitude;

figure;
quiver(u_grid, v_grid, du_norm, dv_norm, 0.5, 'k'); % 0.5 scales arrow length
hold on;

plot(0, 0, 'ro', 'MarkerSize', 5, 'DisplayName', 'Saddle @ (0,0)');
plot(k, 0, 'go', 'MarkerSize', 5, 'DisplayName', 'Stable @ (k,0)');
%plot(3.8891, 1.5844, 'bo', 'MarkerSize', 5, 'DisplayName', 'Saddle @ (u_1^*, v_1^*)');
%plot(1.3661, 2.7565, 'co', 'MarkerSize', 5, 'DisplayName', 'Asymp.Stable @ (u_2^*, v_2^*)')

xlabel('u','FontSize',14);
ylabel('v','FontSize',14);
legend('Location','northeast');
axis([-1 6 -1 6]);
grid on;

%%
A = [0 1; -2 -3];
eig_vals = eig(A);             
lyapunov_exponents = real(eig_vals);  
disp(lyapunov_exponents);

syms t
Phi = expm(A*t);
disp(Phi);


