clear all; clc;


syms lambda a

% a=-5;
b=14;

eqns = lambda^3+(a/6+1)*lambda^2+(-5*a/6+b)*lambda+a*b/6 == 0;

A = 1;
B = (a/6+1);
C = (-5*a/6+b);
D = a*b/6;

Delta = 18*A*B*C*D - 4*B^3*D + B^2*C^2 - 4*A*C^3 - 27*A^2*D^2;

S = solve(eqns, lambda, Maxdegree=3)


%%
a_vals = linspace(6, 14, 300);
b = 14;
hold on
r_list = zeros(200,3);
i=1;

for a = a_vals
    coeffs = [1, (a/6 + 1), (-5*a/6 + b), a*b/6];
    r = roots(coeffs);
    [~, idx] = sort(imag(r));
    r = r(idx);
    r_list(i,:) = r;
    i = i + 1;
    if a == a_vals(1)
%         plot(real(r(1)), imag(r(1)), 'bx', 'MarkerSize', 10);
%         plot(real(r(2)), imag(r(2)), 'rx', 'MarkerSize', 10);
%         plot(real(r(3)), imag(r(3)), 'gx', 'MarkerSize', 10);
%     else
%         plot(real(r(1)), imag(r(1)), 'b.');
%         plot(real(r(2)), imag(r(2)), 'r.');
%         plot(real(r(3)), imag(r(3)), 'g.');
    end
end
% yline(0, 'k');
% xline(0, 'k');
% xlabel('Real Part')
% ylabel('Imaginary Part')
% %title('Root movement as a varies (b = 1)')
% grid on
% axis equal

%%
a = 2;
b_vals = linspace(-5, 5, 100);
hold on

for b = b_vals
    coeffs = [1, (a/6 + 1), (-5*a/6 + b), a*b/6];
    r = roots(coeffs);
    if b == b_vals(1)
        plot(real(r(1)), imag(r(1)), 'x');
        plot(real(r(2)), imag(r(2)), 'x');
        plot(real(r(3)), imag(r(3)), 'x');
    else
        plot(real(r(1)), imag(r(1)), 'b.');
        plot(real(r(2)), imag(r(2)), 'r.');
        plot(real(r(3)), imag(r(3)), 'g.');
    end
end

xlabel('Real Part')
ylabel('Imaginary Part')
title('Root movement as a varies (a = 2)')
grid on
axis equal


%%
[a, b] = meshgrid(-10:0.1:30, -10:0.1:30);

% Coefficients of the characteristic polynomial
a2 = a/6 + 1;
a1 = (-5*a)/6 + b;
a0 = a .* b / 6;

% Compute Routh-Hurwitz condition from row 2 of the Routh table
R1 = ((a2 .* a1) - a0); %./ a2;

% Logical mask for where all Routh-Hurwitz conditions are satisfied
stability = (a2 > 0) & (R1 > 0) & (a0 > 0);

Plot stability region
figure;
imagesc([-10, 20], [-10, 20], stability);  % Plot as image
axis xy;
colormap([1 1 1; 0 0.5 0]);  % white for unstable, green for stable
xlabel('a');
ylabel('b');
title('Stability Region in (a, b) Space');
colorbar('Ticks', [0, 1], 'TickLabels', {'Saddle Point/Unstable', 'Stable'});
grid on;


%%
[a, b] = meshgrid(-10:0.1:30, -10:0.1:30); 

% Coefficients of the characteristic polynomial
a2 = a/6 + 1;
a1 = (-5*a)/6 + b;
a0 = a .* b / 6;

% Discriminant (use element-wise operators!)
Delta = 18 .* a2 .* a1 .* a0 - 4 .* (a2.^3) .* a0 + (a2.^2) .* (a1.^2) - 4 .* (a1.^3) - 27 .* (a0.^2);

% Region with complex conjugate roots
complex_root = Delta < 0;

% Plot
figure;
%figure('Units','normalized','Position',[0.1 0.1 0.6 0.6]);
imagesc([-10,30],[-10,30], complex_root);
axis xy;
colormap([1 1 1; 0 0.5 1]);  % white for real roots, blue for complex roots
xlabel('a');
ylabel('b');
%title('Region with Complex Conjugate Roots in (a, b) Space');
colorbar('Ticks', [0, 1], 'TickLabels', {'All Real', 'Complex Conjugates'});
%legend({'All Real', 'Complex Conjugates'}, 'Location', 'eastoutside');
grid on;
%%
clear all; clc;

%a = linspace(6, 14, step);
k = [2, 14]; %a, b
tspan = linspace(0,500,50000);
%x0 = [-100; 0; -100];
%x3 = [30;30; 30];
x1 = [1e-7;0;0];
x2 = [-1e-7;0;0];

options = odeset('RelTol', 1e-12, 'AbsTol', 1e-12*ones(1,3));
%[t,x] = ode45(@(t,x)asystem(t,x,k), tspan, x0);
[t_1,x_1] = ode45(@(t,x)asystem(t,x,k), tspan, x1, options);
[t_2,x_2] = ode45(@(t,x)asystem(t,x,k), tspan, x2, options);
%[t_3,x_3] = ode45(@(t,x)asystem(t,x,k), tspan, x3);

figure;
hold on;
%plot(t_3, x_3);
% Reset color order to MATLAB default
set(gca, 'ColorOrder', get(groot, 'defaultAxesColorOrder'));

colors = get(gca, 'ColorOrder');

%Plot trajectories (colors are assigned in order: default blue, red, yellow, etc.)
%plot3(x(:,1), x(:,2), x(:,3), 'Color', colors(3,:));
%plot3(x_3(:,1), x_3(:,2), x_3(:,3), 'Color', colors(4,:));
plot3(x_1(:,1), x_1(:,2), x_1(:,3), 'Color', colors(6,:));
plot3(x_2(:,1), x_2(:,2), x_2(:,3), 'Color', colors(1,:));




%plot3(x0(1), x0(2), x0(3), 'ko', 'MarkerSize', 6, 'MarkerFaceColor', colors(3,:));
%plot3(x3(1), x3(2), x3(3), 'ko', 'MarkerSize', 6, 'MarkerFaceColor', colors(4,:));
plot3(x1(1), x1(2), x1(3), 'ko', 'MarkerSize', 6, 'MarkerFaceColor', colors(6,:));
plot3(x2(1), x2(2), x2(3), 'ko', 'MarkerSize', 6, 'MarkerFaceColor', colors(1,:));

xlabel('x');
ylabel('y');
zlabel('z');
axis tight;
legend('(1e-5, 1e-5, 1e-5)','(1e-5, 1e-5, 1e-5)'); %'(1e-5, 1e-5, 1e-5)','(0.1, 0.1, 0.1)', '(1,1,1)',
view(3);


%%
clear all; clc;

% Parameters
b = 14;                     % Fixed b
a_values = 6:0.05:14;       % a range (adjust step size as needed)
tspan_transient = [0 1000]; % Discard transient
tspan_data = [0 500];       % Data collection for maxima
x0 = [1; 1; 1];            % Initial condition

figure; hold on;
xlabel('Parameter $a$', 'Interpreter', 'latex');
ylabel('Local maxima of $x$', 'Interpreter', 'latex');
title('Bifurcation Diagram ($b=14$)', 'Interpreter', 'latex');

% Loop over a values
for i = 1:length(a_values)
    a = a_values(i);
    k = [a, b];
    
    % Discard transient
    [~, x_trans] = ode45(@(t, x) asystem(t, x, k), tspan_transient, x0);
    x0 = x_trans(end, :)'; % Update initial condition
    
    % Collect data for maxima detection
    [t, x_data] = ode45(@(t, x) asystem(t, x, k), tspan_data, x0);
    x_traj = x_data(:, 1); % Extract x-component
    
    % Find local maxima (peaks) in x(t)
    [pks, ~] = findpeaks(x_traj); 
    
    % Plot maxima vs. current a
    if ~isempty(pks)
        plot(a * ones(size(pks)), pks, '.', 'Color', [0.5 0 0], 'MarkerSize', 1);
    end
    
    %fprintf('a = %.2f done\n', a);
end


%%

% Parameters and time vector
a = 6.7; b = 14; k = [a, b];
tspan = linspace(0, 200, 10000);
x0 = [1; 1; 1];  % Initial condition

% Compute LLE
lambda = lyapunovExponent(x0, tspan, k);
fprintf('LLE for a = %.1f, b = %.1f: %.4f\n', a, b, lambda);


