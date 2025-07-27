function dy = tangentSystem(~, y, k)
    a = k(1); b = k(2);
    
    x = y(1:3);     % x = [x; y; z]
    v = y(4:6);     % perturbation vector
    
    dx = zeros(3,1);
    dx(1) = a * (x(2) - (1/16)*x(1)^3 - (1/6)*x(1));
    dx(2) = x(1) - x(2) + x(3);
    dx(3) = -b * x(2);

    J = [ -a*((3/16)*x(1)^2 + 1/6),   a,     0;
           1,                       -1,     1;
           0,                      -b,     0 ];

    dv = J * v;

    dy = [dx; dv];
end
