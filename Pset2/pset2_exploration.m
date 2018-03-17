
% load the data
% times are in 10ths of ms, 0.0001s
% only work with neurons [1, 4, 11, 15, 26, 51, 80, 84, 96, 105]
load('retinaData.mat');

selected_neurons = {};

indices = [1, 4, 11, 15, 26, 51, 80, 84, 96, 105];

for i=1:length(indices)
    selected_neurons{i} = retinaData.spikes{indices(i)};
end

%% raster

figure; hold on

% get the first minute
frame = 60 / .0001;
for i=1:length(selected_neurons)
    neur = selected_neurons{i};
    neur = neur(neur <= frame);
    for iid = 1:length(neur)
        spkx=[neur(iid) neur(iid)];
        spky = [0 1] + i;
        line(spkx,spky,'LineWidth',.5);
    end
end


%% strf 2a

% get strf for neuron 1
stimulus = retinaData.stimulusFrames;

times = retinaData.stimulusFrameTimes;

%strf_by_time = zeros(40, 40, 1);

% create an array of lag times
lags = 300:-10:-10;
lags = lags / .1;
% cite this
figure('rend','painters','pos',[10 10 900 600]); hold on

for i=1:length(lags)
    m = getSTRF(selected_neurons{1}, stimulus, times, lags(i));
    subplot(6,6,i);
    pcolor(m);
    colormap('gray')
    title(sprintf('Lag = %d ', lags(1) / 10))
end


suptitle('STRFs for Neuron 1 by Time Lag')

%% 2b, get STRFs for all neurons, plot and save

strf_by_neuron = {};

for i=1:length(selected_neurons)
    f = figure('rend','painters','pos',[10 10 900 600]); hold on
    set(f, 'visible', 'off');
    strf_by_neuron{i} = zeros(40,40,1);
    for j=1:length(lags)
        % repeating some code because it's small enough that it doesnt make
        % sense to make a new function
        m = getSTRF(selected_neurons{i}, stimulus, times, lags(j));
        strf_by_neuron{i}(:,:,j) = m;
        subplot(6,6,j);
        pcolor(m);
        colormap('gray')
        title(sprintf('Lag = %d ', lags(1) / 10))
    end
    suptitle(sprintf('Neuron %d STRF', i));
    fig_name = sprintf('Neuron_%d STRF', i);
    saveas(f, fig_name, 'fig');
end





%% 3

% use one value of sigma
% amplitude is constant multiple 
% A + Bg(x,y) so that you have an intercept
% provide initial guess that is close to the center 
% fit using the most pronounced strf with no parameters set 
% initial guess for constant should be .5
% x is meshgrid(1:40)
% y is the strf values
% reshape meshgrid and strf values so that they match


% creating meshgrid and reshaping 
[a1, a2] = meshgrid(1:40);
r(:,1) = reshape(a1(:,:), [40^2 1]);
r(:,2) = reshape(a2(:,:), [40^2 1]);

params = {};

for i=1:length(strf_by_neuron)
    params{i} = BuildModel(strf_by_neuron{i}, r);
end



%% Plotting


lags_in_ms = lags / 10;

ind_80 = find(lags_in_ms == 80)

for i=1:length(strf_by_neuron)
    figure; hold on
    p = params{i};
    v(1) = p{1};
    v(2) = p{2};
    v(3) = p{3};
    v(4) = p{4}(ind_80);
    v(5) = p{5}(ind_80);
    preds = Gauss2D(v, r);
    strf = strf_by_neuron{i};
    strf = strf(:,:,ind_80);
    subplot(3,1, 1)
    pcolor(strf);
    subplot(3,1,2)
    surf(reshape(r(:,1), 40, 40), reshape(r(:,2), 40, 40), reshape(preds, 40,40));
    subplot(3,1,3)
    plot(lags_in_ms, p{4})
    ylabel('Amplitudes')
    xlabel('Lag Time')
    suptitle(sprintf('Neuron %d', i))

end


%% Final Questions 


fprintf('The model did not seem to accurately fit every neuron, at least')
fprintf('At the 80ms lag STFR. For some it displayed a peak where there was')
fprintf('none in the STRF. On the other hand, it did seem to fit well for others')


widths = [];

for i = 1:length(params)
    widths = [widths params{i}{2}];
end

%% Average Width

fprintf('The average width was %2.4f', mean(widths))

%% ON and OFF

count = 0;
for i=1:length(strf_by_neuron)
    [best, best_strf, best_ind] = GuessCenter(strf_by_neuron{i});
    % code is repeated only because it doesn't make sense to make a
    % separate function just for selecting data from parameters
    p = params{i};
    v(1) = p{1};
    v(2) = p{2};
    v(3) = p{3};
    strf = strf_by_neuron{i};
    strf = strf(:,:,best_strf);
    strf = strf(:);
    mean_over_whole = mean(strf);
    rf_vals = []
    for j = 1:length(r(:,1))
        % find coordinates within RF to get mean over RF
        coord = r(i,:);
        dist = sqrt(sum((coord - v(1:2)).^2));
        if dist <= (v(3) / 2)
            rf_vals = [rf_vals strf(i)];
        end
    end
    mean_over_rf = mean(rf_vals);
    if mean_over_rf <= mean_over_whole
        count = count + 1;
    end
end

    
% im getting a lot of nan values for mean but I don't have time to figure
% out why so im just going to go ahead

% 10, 9, 8, 7, 4, 3 are OFF from the plots

fprintf('The percentage of neurons that are OFF is %d ', 60)

offs = [10 9 8 7 4 3];

ons = [6 5 2 1];

off_avg = GetAvg(offs, lags_in_ms, params);
on_avg = GetAvg(ons, lags_in_ms, params);
fprintf('The average time until peak for OFF neurons was %3.2f ms', off_avg);
fprintf('The average time until peak for ON neurons was %3.2f ms', on_avg);





        