function trial_by_direction = SortTrialsByDirection(directions, neuron)

% create a map to associate directions to trials




trial_by_direction = containers.Map(double(5), [0 0 0]);
trial_by_direction.remove(5);
for i=1:length(neuron)
    trial = neuron{i};
    dir = directions(i);
    % check if a direction has already been added to the map
    if ~trial_by_direction.isKey(dir)
        newcell = {};
        newcell{1} = trial;
        trial_by_direction(dir) = newcell;
    else
        cell = trial_by_direction(dir);
        ind = length(cell) + 1;
        cell{ind} = trial;
        trial_by_direction(dir) = cell;
    end
    
end

end