function [spike_times, V, t, c] = model(A)

% set timestep to 1 ms
timestep = .001;

t = 0:timestep:.5;
[c, t] = Current(t, A);

V = [];
spike_times = [];

% deal with the initial conditions
V(1) = 0;
dv_1 = LIF(c(1), V(1));

[new_V, has_spiked] = update(dv_1, V(1), timestep, false);

V(2) = new_V;

if has_spiked
    spike_times = [spike_times t(1)]
end

% run the model for the rest of the times
for i=2:length(c)
    
    curr = c(i);
    dv = LIF(curr, V(i));
    [new_V, has_spiked] = update(dv, V(i), timestep, has_spiked);
    V(i+1) = new_V;
    if has_spiked
        spike_times = [spike_times t(i)];
    end
end

% fill in initial conditions for t and c so that they can be plotted
t = [0 t];
c = [0 c];


end