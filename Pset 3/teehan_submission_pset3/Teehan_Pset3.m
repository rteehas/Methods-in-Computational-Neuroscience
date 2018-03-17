%% 1.1



[spike_times, V, t, current] = model(3);

figure; hold on

subplot(2,1,1)
plot(t, V)
xlabel("Time (s)")
ylabel("Membrane Voltage")

subplot(2,1,2)
plot(t,current)
ylabel("Current")
xlabel("Time (s)")

%% 1.1 Question

% period is 1/freq, number of cycles is time/period
num_cycles = .5/(1/40);
firing_rate = length(spike_times)/.5;
fprintf("The firing rate was: %2.1f spikes/sec\n", firing_rate);
fprintf("It fires %1.3f spikes per cycle\n", length(spike_times) / num_cycles);

%% 1.2 Johnson Plot

amplitudes = linspace(0, 6, 100);

times = {};
rate_per_cycle = [];

for i=1:length(amplitudes)
    A = amplitudes(i);
    [spike_times, V, t, c] = model(A);
    total = length(spike_times);
    % or total / (1/40)
    cycle_rate = total / num_cycles;
    rate_per_cycle = [rate_per_cycle cycle_rate];
    times{i} = spike_times;
end

figure; hold on
plot(amplitudes, rate_per_cycle)
ylabel("Impulses per cycle")
xlabel("Amplitude")

%% Raster

figure; hold on

diff = 0.0600;

for i=1:length(times)
    trial = times{i};
    
    for iid = 1:length(trial)
        hold on
        spkx=[trial(iid) trial(iid)];
        spky = [0 diff] + i*diff;
        line(spkx,spky,'LineWidth',1);
    end
end
hold on
xlabel("Time (s)")
ylabel("Amplitudes")
ylim([0 6.5])

%% Questions


fprintf("The model is successful insofar as it is able model a simplified ")
fprintf("version of action potentials, in which a spike occurs at a certain ")
fprintf("threshold and resets to a lower value.\n\n")
fprintf("The major failure in this model, however, is that it does not deal ")
fprintf("with adaptation since there is no memory of past spikes.\n")
fprintf("Another way that it differs from natural behavior, and also a feature ")
fprintf("of real neurons that it does not capture is that real neurons get ")
fprintf("inputs from a number of different neurons in the network.")
fprintf(" A single sinusoidal input current may not be able to properly model that.")
fprintf(" Looking at the raster plot, at a certain point the spikes happen")
fprintf(" almost entirely at the same set of times, which may indicate a kind")
fprintf(" of stable regular firing time.")






    