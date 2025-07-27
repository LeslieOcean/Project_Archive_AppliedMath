function lambda = lyapunovExponent(x0, tspan, params)
    d0 = 1e-6;             
    dt = tspan(2) - tspan(1);
    N = length(tspan);
    
    y0 = [x0; d0*randn(3,1)];           % [state; perturbation vector]
    y0(4:6) = y0(4:6) / norm(y0(4:6)) * d0;  %normalize
    
    Y = zeros(N, 6);
    Y(1,:) = y0';
    sum_log = 0;
    steps = 0;
    
    for i = 2:N
        [~, y] = ode45(@(t,y) tangentSystem(t, y, params), [tspan(i-1), tspan(i)], Y(i-1,:)');
        y_end = y(end,:)';
        
        x = y_end(1:3);
        v = y_end(4:6);
        d = norm(v);
        
        if d > 0
            sum_log = sum_log + log(d / d0);
            steps = steps + 1;
        end
        
        v = v / d * d0;
        Y(i,:) = [x; v]';
    end
    
    lambda = sum_log / (steps * dt);
end