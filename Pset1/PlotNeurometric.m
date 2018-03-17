
% NOTE initialize a figure before running this

function threshold = PlotNeurometric(rates_by_direction, preferred, num_samples)

keys = rates_by_direction.keys;
bootstrap = zeros(length(keys),1);
ang_diff = zeros(length(keys),1);

for i=1:length(keys)
    k = keys{i};
    key_firing = rates_by_direction(k);
    pref_firing = rates_by_direction(preferred);
    key_sample = randsample(key_firing, num_samples, true);
    preferred_sample = randsample(pref_firing, num_samples, true);
    bootstrap(i) = sum(key_sample <= preferred_sample) / num_samples;
    ang_diff(i) = abs(preferred - k);
end

scatter(ang_diff, bootstrap, [], 'black', 'filled');
xlabel('Angular difference (°)')
ylabel('Probability')

% to find the threshold we find the value for bootstrap closest to .75
% in the beginning of the bootstrap values
% ideally I would interpolate or find the intertersection of the graph with
% the line at y = .75, but I can't figure out how to do that

[val, id] = min(abs(bootstrap - .75));

threshold = ang_diff(id);
end
