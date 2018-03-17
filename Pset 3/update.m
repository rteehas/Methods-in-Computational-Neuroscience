% spiked is a boolean 
function [new_V, has_spiked] = update(dV, V, timestep, spiked)


spike_v = 4;
thresh_v = 1;
reset_v = -1;

has_spiked = false;

new_V = V + dV * timestep;

% check if the neuron has just spiked and reset
if spiked
    new_V = reset_v;
    has_spiked = false;

% if not, check if the voltage is greater than the threshold
elseif new_V > thresh_v
    new_V = spike_v;
    has_spiked = true;
end
end

