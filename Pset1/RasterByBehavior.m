function RasterByBehavior(neuron, behaviors, behavior)
for i=1:length(neuron)
    if behaviors(i) == behavior
        trial = neuron{i};
        for iid=1:length(trial);
            spkx=[trial(iid) trial(iid)];
            spky = [0 1] + i;
            line(spkx,spky,'LineWidth',1);
        end
    end

end
% add lines to separate pre-trial, sample, delay and test periods
line([0 0], get(gca, 'ylim'), 'color', 'black');
line([.65 .65], get(gca, 'ylim'), 'color', 'black');
line([1.65 1.65],get(gca, 'ylim'), 'color', 'black');

xlabel('Time (sec')
end
        