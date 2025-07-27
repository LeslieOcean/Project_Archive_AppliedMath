function dphi = kuramoto_3oscillator_1(t, phi, domega, k)

dphi = [
    domega(1) + k/3*(-2*sin(phi(1))-sin(phi(1)+phi(2))+sin(phi(2)));
    domega(2) + k/3*(sin(phi(1))-sin(phi(1)+phi(2))-2*sin(phi(2)));
];