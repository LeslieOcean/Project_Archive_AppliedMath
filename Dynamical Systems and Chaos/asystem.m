function dx = asystem(t, x, k)

dx = [
    k(1)*(x(2)-1/16*x(1)^3-1/6*x(1));
    x(1)-x(2)+x(3);
    -k(2)*x(2);
];