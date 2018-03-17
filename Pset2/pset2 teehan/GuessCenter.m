% gives a rough guess of the center of the gaussian 

function [best, best_strf, best_ind] = GuessCenter(neuron_strf)

best = 0;

for i=1:32
    strf = neuron_strf(:,:,i);
    diffs = strf - mean(strf(:));
    [d, ix] = min(diffs(:));
    if d < best
        best = d;
        best_ind = ix;
        best_strf = i;
    end
end

end
