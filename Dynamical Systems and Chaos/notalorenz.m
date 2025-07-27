function dx = notalorenz(t, x, beta)

dx = [
    beta(1)*(x(2)-x(1));
    beta(2)*x(1)-x(2)+x(1)*x(3);
    x(1)*x(2)-beta(3)*x(3);
];