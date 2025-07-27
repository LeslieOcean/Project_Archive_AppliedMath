function dphi = kuramoto(t, theta, omega, k)

dphi = [
    omega(1) + k/2*sin(theta(2)-theta(1));
    omega(2) + k/2*sin(theta(1)-theta(2));
];