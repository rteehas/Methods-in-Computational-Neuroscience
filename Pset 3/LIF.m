function dV = LIF(current, V)


rest_V = 0;

% 15ms
tau = 0.015;

dV = (-(V - rest_V) + current) / tau;

end