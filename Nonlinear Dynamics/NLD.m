clear all; clc;

omega = [1; 0.5]; %w1, w2
k = 0.5;
fp = asin((omega(1)-omega(2))/k);
theta0 = [fp-0.1; fp]; %ICs
tspan = linspace(0,100,5000);

%options = odeset('RelTol', 1e-12, 'AbsTol', 1e-12*ones(1,2));
[t,theta] = ode45(@(t,theta)kuramoto(t,theta,omega,k), tspan, theta0);
phi = theta(:,1) - theta(:,2);

figure(1);
plot(t,sin(theta(:,1)));
hold on;
plot(t,sin(theta(:,2)));
legend('sin(\theta_1)','sin(\theta_2)')

figure(2);
plot(t,cos(phi));
legend('cos(\phi)')

figure(3);
plot(t,phi);
legend('\phi');

%%
clear all; clc;

itera =100;

omega = [1; 0.5]; %w1, w2
%k = linspace(0, 2, itera);
k = 0.3;
IC = 2*pi*rand(itera,2);
phi_t_list = zeros(itera);

% for j=1:1:itera
%     phi_t = 0;
for i=1:1:itera
    %fp = asin((omega(1)-omega(2))/k(i));
    
    theta0 = [IC(i,1); IC(i,2)]; %ICs
    tspan = linspace(0,100,5000);

    %options = odeset('RelTol', 1e-12, 'AbsTol', 1e-12*ones(1,2));
    [t,theta] = ode45(@(t,theta)kuramoto(t,theta,omega,k), tspan, theta0);
    phi = theta(:,1) - theta(:,2);
    %phi_t = phi_t + cos(phi(end));
    figure(1);
    scatter(i, cos(phi(end)));
    hold on
end
% phi_t = phi_t/itera;
% phi_t_list(j) = phi_t;
% end

% figure();
% plot(k, phi_t_list);

%%
clear all; clc;

k = 0.25;
itera = 40;

omega1 = linspace(0, 1, itera);
omega2 = 0.5;
domega = omega1 - omega2;
IC = 2*pi*rand(itera,2);
phi_dot = zeros(itera,1);

for i = 1:1:itera %for different domega
    phi_t = 0;
    for j = 1:1:itera %for different ICs
        omega = [omega1(i), omega2];
        theta0 = [IC(j,1); IC(j,2)]; %ICs
        tspan = linspace(0,100,5000);
        
        %options = odeset('RelTol', 1e-12, 'AbsTol', 1e-12*ones(1,2));
        [t,theta] = ode45(@(t,theta)kuramoto(t,theta,omega,k), tspan, theta0);
        phi = theta(:,1) - theta(:,2);
        phi_t = phi_t + domega(i) - k*sin(phi(end)); 
    end
    phi_t = phi_t/itera;
    phi_dot(i) = phi_t;
end

figure();
plot(domega, phi_dot)

%%
clear all; clc;


itera = 10;
k = linspace(0, 0.6, itera);

omega1 = linspace(0, 1, itera);
omega2 = 0.5;
domega = omega1 - omega2;
IC = 2*pi*rand(itera,2);
phi_dot = zeros(itera,itera);

for m = 1:1:itera %for different k
    for i = 1:1:itera %for different domega
        phi_t = 0;
        for j = 1:1:itera %for different ICs
            omega = [omega1(i), omega2];
            theta0 = [IC(j,1); IC(j,2)]; %ICs
            tspan = linspace(0,100,5000);
            
            %options = odeset('RelTol', 1e-12, 'AbsTol', 1e-12*ones(1,2));
            [t,theta] = ode45(@(t,theta)kuramoto(t,theta,omega,k(m)), tspan, theta0);
            phi = theta(:,1) - theta(:,2);
            phi_t = phi_t + domega(i) - k(m)*sin(phi(end)); 
        end
        phi_t = phi_t/itera;
        phi_dot(m,i) = phi_t;
    end
end


figure();
h = heatmap(domega, k, phi_dot);







