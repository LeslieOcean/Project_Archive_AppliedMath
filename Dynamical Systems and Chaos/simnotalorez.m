clear all; clc;

%a = linspace(6, 14, step);
beta = [10, 28, 8/3]; %a, b
tspan = linspace(0,10,5000);
x0 = [1; 1; 1];
x1 = [0.5;0.5;0.5];
x2 = [1e-5;1e-5;1e-5];
x3 = [5; 5; 5];

%options = odeset('RelTol', 1e-12, 'AbsTol', 1e-12*ones(1,3));
[t,x] = ode45(@(t,x)notalorenz(t,x,beta), tspan, x0);
% [t_1,x_1] = ode45(@(t,x)asystem(t,x,k), tspan, x1);
% [t_2,x_2] = ode45(@(t,x)asystem(t,x,k), tspan, x2);
% [t_3,x_3] = ode45(@(t,x)asystem(t,x,k), tspan, x3);

figure;
hold on;
%plot(t_3, x_3);
% Reset color order to MATLAB default
set(gca, 'ColorOrder', get(groot, 'defaultAxesColorOrder'));

colors = get(gca, 'ColorOrder');

%Plot trajectories (colors are assigned in order: default blue, red, yellow, etc.)
plot3(x(:,1), x(:,2), x(:,3), 'Color', colors(3,:));
%plot3(x_3(:,1), x_3(:,2), x_3(:,3), 'Color', colors(4,:));
%plot3(x_1(:,1), x_1(:,2), x_1(:,3), 'Color', colors(6,:));
% plot3(x_2(:,1), x_2(:,2), x_2(:,3), 'Color', colors(1,:));




plot3(x0(1), x0(2), x0(3), 'ko', 'MarkerSize', 6, 'MarkerFaceColor', colors(3,:));
% plot3(x3(1), x3(2), x3(3), 'ko', 'MarkerSize', 6, 'MarkerFaceColor', colors(4,:));
% %plot3(x1(1), x1(2), x1(3), 'ko', 'MarkerSize', 6, 'MarkerFaceColor', colors(6,:));
% plot3(x2(1), x2(2), x2(3), 'ko', 'MarkerSize', 6, 'MarkerFaceColor', colors(1,:));

xlabel('x');
ylabel('y');
zlabel('z');
axis tight;
%legend('(5,5,5)','(1e-5, 1e-5, 1e-5)'); %'(1e-5, 1e-5, 1e-5)','(0.1, 0.1, 0.1)', '(1,1,1)',
view(3);


%%

syms lambda

sigma = 10;
r = 28;
b=8/3;

%eqns = lambda^3 - (-sigma-1-b)*lambda^2 + (-270+88/3)*lambda - 1960/8;
eqns = lambda^2+11*lambda-270;

S = solve(eqns==0, lambda, Real=true);

%%
[x, y, z] = meshgrid(linspace(0, 200, 100), ...
                     linspace(0, 200, 100), ...
                     linspace(0, 200, 100));

V = 28*x.^2 + 10*y.^2 + 10*(z - 56).^2;
V_dot = 28*x.^2 + y.^2 + (8/3)*(z - 28).^2 - (8/3)*28^2- 2*x.*y.*z;

level = 1e4;

figure;
p = patch(isosurface(x, y, z, V_dot, level));
isonormals(x, y, z, V, p)
set(p, 'FaceColor', 'red', 'EdgeColor', 'none');

camlight; lighting gouraud;
xlabel('x'); ylabel('y'); zlabel('z');
title('Isosurface of V(x,y,z) = 0');
axis equal;
grid on;
view(3);

%%
A = [-10, 10, 0;
      28, -1, 0;
      0,  0, -8/3];

eigvals = diag(D);
eigvecs = V;

unstable_idx = find(real(eigvals) > 0);
unstable_vec = eigvecs(:, unstable_idx);

stable_idx = find(real(eigvals) < 0);
stable_vecs = eigvecs(:, stable_idx);

