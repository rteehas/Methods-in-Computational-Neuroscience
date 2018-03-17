
load('S1_Ideal_Observer_Analysis.mat');


%% 3.1a


% preparing the data
amplitudes = Data.stimuli{1}(:,2,:);
directions = Data.stimuli{1}(:,3,:);

new_directions = [];
neuron_1 = {};
neuron_2 = {};
for i = 1:length(Data.spikes{1})
    if amplitudes(i) == 700
        neuron_1{length(neuron_1) + 1} = Data.spikes{1}{i};
        neuron_2{length(neuron_2) + 1} = Data.spikes{2}{i};
        new_directions = [new_directions directions(i)];
    end
end

% make a colormap
colormap = GetColormap(directions);



figure; hold on


subplot(1,2,1)
RasterByDirection(new_directions, neuron_1, colormap)
title('Neuron 1');


subplot(1,2,2)
RasterByDirection(new_directions, neuron_2, colormap)
title('Neuron 2');



%% 3.1b


figure; hold on

% sort trials by direction

trial_by_direction_1 = SortTrialsByDirection(new_directions, neuron_1);
trial_by_direction_2 = SortTrialsByDirection(new_directions, neuron_2);

% now get the avg firing rates for each direction

% now plot and get the preferred direction and firing rates sorted
% by direction
subplot(2,1,1)
[preferred_1, rates_by_direction_1] = ...
    GetFiringRatesAndAvg(trial_by_direction_1);
title('Neuron 1')

subplot(2,1,2)
[preferred_2, rates_by_direction_2] = ...
    GetFiringRatesAndAvg(trial_by_direction_2);
title('Neuron 2')



%% 3.1c
figure; hold on

subplot(2,1,1)
threshold_1 = PlotNeurometric(rates_by_direction_1, preferred_1, 1000)
title('Neuron 1');

subplot(2,1,2)
threshold_2 = PlotNeurometric(rates_by_direction_2, preferred_2, 1000)
title('Neuron 2')


fprintf('The threshold for both is around 22.5 \n')
fprintf('The second neuron is abnormal because it dips significantly \n')
fprintf(' below the threshold after it reaches it')

%% 3.1d

[preferred_1, new_rates_1] = ResortByDirection(rates_by_direction_1);
[preferred_2, new_rates_2] = ResortByDirection(rates_by_direction_2);

figure; hold on

subplot(2,1,1)
threshold_1 = PlotNeurometric(new_rates_1, preferred_1, 10000)

subplot(2,1,2)
threshold_2 = PlotNeurometric(new_rates_2, preferred_2, 10000)

fprintf('The thresholds are: \n')
fprintf('Neuron 1: %2.1f \n', threshold_1)
fprintf('Neuron 2: %2.4f \n', threshold_2)


%% 3.2

load('choiceData.mat')


%% 3.2a

behaviors_1 = choiceData.behavioralReport{1};
behaviors_2 = choiceData.behavioralReport{2};
neuron_1 = choiceData.spikes{1};
neuron_2 = choiceData.spikes{2};

figure;hold on

RasterByBehavior(neuron_1, behaviors_1, 1)

figure; hold on
RasterByBehavior(neuron_1, behaviors_1, 2)

figure; hold on 
RasterByBehavior(neuron_2, behaviors_2, 1)

figure; hold on
RasterByBehavior(neuron_2, behaviors_2, 2)

%% 3.2b

p_1 = ChoiceProbabilities(neuron_1, behaviors_1, 1000);
p_2 = ChoiceProbabilities(neuron_2, behaviors_2, 1000);


%% 3.2b
fprintf('The choice probabilities are: \n')
fprintf('Neuron 1 \n')
fprintf('Pretrial: %2.4f \n', p_1(1)) 
fprintf('Sample: %2.4f \n', p_1(2))
fprintf('Delay: %2.4f \n', p_1(3))
fprintf('Test: %2.4f', p_1(4))
fprintf('\n \n')
fprintf('Neuron 2 \n')
fprintf('Pretrial: %2.4f \n', p_2(1)) 
fprintf('Sample: %2.4f \n', p_2(2))
fprintf('Delay: %2.4f \n', p_2(3))
fprintf('Test: %2.4f \n', p_2(4))

fprintf('The choice probabilities in the pre-trial intervals are')
fprintf(' %2.4f and %2.4f \n', [p_1(1) p_2(1)])
fprintf('This means that even among the ambiguous trials, the ones which')
fprintf(' were coded with a 1 had a different distribution than those\n')
fprintf('that were coded with a 2\n')
fprintf('In other words, the codings for ambiguous data were not randomly')
fprintf('assigned')






