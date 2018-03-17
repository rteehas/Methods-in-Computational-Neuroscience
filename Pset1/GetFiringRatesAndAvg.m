
% NOTE, initialize a figure before this

function [preferred, rates_by_direction] = ...
    GetFiringRatesAndAvg(trial_by_direction)


x = zeros(16, 1);
y = zeros(16,1);
err = zeros(16, 1);
keys = trial_by_direction.keys;

% get the rates by direction to be used later
rates_by_direction = containers.Map(double(5), [0 0 0]);
rates_by_direction.remove(5);
for i=1:length(keys)
    k = keys{i};
    trial = trial_by_direction(k);
    rates = zeros(length(trial), 1);
    for j=1:length(trial)
        t = trial{j};
        if ~isempty(t);
            rates(j) = length(t) / t(length(t));
        end 
    end
    rates_by_direction(k) = rates;
    x(i) = k;
    y(i) = mean(rates);
    % get the standard error of the mean
    err(i) = std(rates) / sqrt(length(rates));
end

% get argmax to get preferred direction 
[argvalue, argmax] = max(y);

preferred = x(argmax);

% plot tuning curve
errorbar(x, y, err)
xlabel('Direction');
ylabel('Mean Firing Rate');

end
