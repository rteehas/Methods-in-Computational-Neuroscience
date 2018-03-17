
load('mtSpikeTimes.mat');
load('S1_Ideal_Observer_Analysis.mat');

%% 

% 1.1
% raster plot code, to be made into function 
figure; hold on
subplot(3,,3)
for i = 1:length(mtSpikeTimes)
    hold on
    trial = mtSpikeTimes{i}
    
    for iid = 1:length(trial)
        hold on
        spkx=[trial(iid) trial(iid)]
        spky = [0 1] + i
        line(spkx,spky,'LineWidth',1)
    end
end
xlabel('Time (sec)')
%% 

%1.2 
% include function that takes bin size as variable and plot a couple of
% them
% make this a function with different bin sizes
% figure;
% 
% edges = linspace(0, .6, 60); 
% psth = zeros(1,length(edges));
% 
% for j = 1:length(mtSpikeTimes)
%     psth = psth + histc(mtSpikeTimes{j}, edges);
% end
% 
% avg = psth / length(mtSpikeTimes)
% plot(edges, avg);
% xlim([0 .6]);
% xlabel('Time (sec)');
% ylabel('# of spikes');

subplot(3,2,1)
PlotPSTH(mtSpikeTimes, 10)

subplot(3,2,2)
PlotPSTH(mtSpikeTimes, 20)


%% 2.1

%mtspike data is time in seconds, divide by ms time to get number 

rates = zeros(184, 1);
counts = zeros(184, 1);
range = .3;
d = [];
for i = 1:length(mtSpikeTimes)
    if length(mtSpikeTimes{i}) > 0
        cutoff = mtSpikeTimes{i}(1) + range;
        count = sum(mtSpikeTimes{i} <= cutoff);
        counts(i) = count;
        rates(i) = count / range;
        vals = mtSpikeTimes{i}(mtSpikeTimes{i} <= cutoff);
        new = diff(vals);
        d = [d new];
    end
end

% get mean, std, and coefficient of variation
mean_val = mean(rates)
std_val = std(rates)
coef_var = std_val / mean_val

% fano factor
fano = var(counts) / mean(counts)

fprintf('The Fano Factor of the counts is: %4.5f \n', fano)



%% 2.2
% isi histogram
close all;
figure; hold on 
%normalize histogram
[nums, cents] = hist(d, 100);
nums = nums / nums(1);
bar(cents, nums);


% fix curve fitting
x = linspace(0, .5, 100);
y = exp((-x)*mean_val);
plot(x,y, 'color', 'red');



%% test functions

PlotPSTH(mtSpikeTimes, 7)




