% NOTE: Times are in 10th of a ms
% should only have one frame of stimulus for each spike
% find frame related ot each spike w lag
% mean(arg, 3)
% get separate STRF for each lag 
% each stimulus frame i is of the form retinaData.stimulusFrames(:,:,i)

function strf = getSTRF(neuron, stimulus, stim_times, lag)

stims = zeros(40,40,1);
for i =1:length(neuron)
    time = neuron(i);
    new_time = time - lag;
    index = find(stim_times == new_time);
    stims(:,:,index) = stimulus(:,:,index);
end

strf = mean(stims, 3);
end