
% Creates a rasterplot from an array of directions, a cell of trials 
% associated with a neuron, and a map that tells you which color to use for 
% each direction. NOTE: a figure needs to be initialized before running 
% the function 

function RasterByDirection(directions, neuron, colormap)



for i=1:length(neuron)
    trial_neuron_1 = neuron{i};
    
    for iid=1:length(trial_neuron_1)
        spkx=[trial_neuron_1(iid) trial_neuron_1(iid)];
        spky = [0 1] + i;
        col = colormap(directions(i));
        line(spkx,spky,'Color', col, 'LineWidth',1);
    end
end
xlabel('Time (sec)');








