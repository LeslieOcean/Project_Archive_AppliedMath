function dphi = kuramoto_3oscillator_2(t, theta, domega, k)

dphi = [
    domega(1) + k/4*(sin(theta(2)-theta(1)));
    domega(2) + k/4*(sin(theta(1)-theta(2))+sin(theta(3)-theta(2)));
    
];