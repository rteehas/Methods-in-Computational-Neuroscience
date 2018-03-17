
% Creates a rasterplot from an array of directions, a cell of trials 
% associated with a neuron, and a map that tells you which color to use for 
% each direction. NOTE: a figure needs to be initialized before running 
% the function 

function RasterByDirection(directions, data, colormap, times)



for i=1:length(data(1,:,1))
    dir = directions(i);
    col = colormap(dir);
    for j=i:length(data(1,1,:))
        trial = data(:,i,j);
    
        for iid=1:length(trial)
            hold on
            if trial(iid) == 1
                t = times(iid);
                spkx=[t t];
                spky = [0 1] + i;
                col = colormap(directions(i));
                line(spkx,spky,'Color', col, 'LineWidth',1);
            end
        end
    end
end
xlabel('Time (sec)');








