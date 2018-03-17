function avg_time = GetAvg(indices, lags_in_ms, params)


times = [];
for i=1:length(indices)
    a = params{i}{4};
    [d, ix] = max(a);
    times = [times lags_in_ms(ix)];
end

avg_time = mean(times);
end