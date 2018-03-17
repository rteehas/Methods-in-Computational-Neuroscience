function probs = ChoiceProbabilities(neuron, behaviors, num_samples)
beh_1_trials = {};
beh_2_trials = {};

for i=1:length(neuron)
    if behaviors(i) == 1
        beh_1_trials{length(beh_1_trials) + 1} = neuron{i};
    else if behaviors(i) == 2
            beh_2_trials{length(beh_2_trials) + 1} = neuron{i};
        end
    end
end

rates_1 = RatesByPeriod(beh_1_trials);
rates_2 = RatesByPeriod(beh_2_trials);

probs = [];
for i=1:length(rates_1)
    sample_1 = randsample(rates_1{i}, num_samples, true);
    sample_2 = randsample(rates_2{i}, num_samples, true);
    prob = sum(sample_1 <= sample_2) / num_samples;
    probs = [probs prob];
end


end
