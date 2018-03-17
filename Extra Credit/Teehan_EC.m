load('spikes.mat')


%% get the values for each spike train

% you can rerun this with larger values of q, but it took too long to run
% for me to be able to get the range shown in the graph in the paper. 

qs = 0:25:500;

f = GetFractions(spikes,qs);
figure;
plot(qs, f);
xlabel('Q values')
ylabel("Accuracy")