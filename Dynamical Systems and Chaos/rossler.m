function dx = rossler(t, x, beta)

dx = [
    -x(2)-x(3);
    x(1)+beta(1)*x(2);
    beta(2)+x(3)*(x(1)-beta(3));
];