function [b, dev, stats, aics, model] = addHistory(model, hist, spikes)

aics = [];

for i=1:20
    model = [model hist(:,i)];
    [b, dev, stats] = glmfit(model, spikes, 'poisson');
    aic = dev + 2*length(model);
    aics = [aics aic];
    
end

end