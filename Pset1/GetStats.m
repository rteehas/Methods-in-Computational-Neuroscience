function [mean_val, std_val, coef_var, fano, d] = GetStats(cell_array, range)
rates = zeros(length(cell_array), 1);
counts = zeros(length(cell_array), 1);
d = [];
for i = 1:length(cell_array)
    if length(cell_array{i}) > 0
        cutoff = cell_array{i}(1) + range;
        count = sum(cell_array{i} <= cutoff);
        counts(i) = count;
        rates(i) = count / range;
        vals = cell_array{i}(cell_array{i} <= cutoff);
        new = diff(vals);
        d = [d new];
    end
end

% get mean, std, and coefficient of variation
mean_val = mean(rates);
std_val = std(rates);
coef_var = std_val / mean_val;

% fano factor
fano = var(counts) / mean(counts);
end
