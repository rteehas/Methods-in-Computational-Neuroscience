function [preferred, new_rates] = ResortByDirection(rates)
keys = rates.keys;

new_rates = containers.Map(double(5), [0 0 0]);
new_rates.remove(5);

for i=1:length(keys)
    dir = keys{i};
    
    % the description in the pset seems to indicate that we take the
    % direction mod 180 and add 90
    
    parallel = dir + 180;
    if rates.isKey(parallel)
        dir_vals = rates(dir);
        par_vals = rates(parallel);
        dir_vals = dir_vals';
        par_vals = par_vals';
        final = [dir_vals par_vals];
        final = final';
        
        new_k = mod(parallel, 180);
        new_rates(new_k) = final;
    end
end

new_keys = new_rates.keys;
firing = [];
dirs = [];
for i=1:length(new_keys)
    k = new_keys{i};
    dir_fire = new_rates(k);
    firing = [firing mean(dir_fire)];
    dirs = [dirs k];
end

[argvalue, argmax] = max(firing);
preferred = dirs(argmax);
end










