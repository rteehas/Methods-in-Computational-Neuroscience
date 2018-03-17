% builds a model using lsqcurvefit

function params = BuildModel(neuron_strf, coords)

% get initial guesses and the strf we should use to fit 
% the center and width

params = {};

[best, best_strf, best_ind] = GuessCenter(neuron_strf);

cent = coords(best_ind,:);

first_test = neuron_strf(:,:,best_strf);

% use some initial guesses for amplitude, sigma, and b
vals = ...
    lsqcurvefit(@Gauss2D, [cent(1), cent(2), 10, 1, .5], coords, first_test(:));

params{1} = vals(1);
params{2} = vals(2);
params{3} = vals(3);

amplitudes = [];
intercepts = [];

for i=1:32
    ydata = neuron_strf(:,:,i);
    ydata = ydata(:);
    final_terms = ...
        lsqcurvefit(@(x0, xdata) Gauss2D([vals(1), vals(2), vals(3), x0(1), x0(2)], xdata), [1, .5], coords, ydata);
    amplitudes = [amplitudes final_terms(1)];
    intercepts = [intercepts final_terms(2)];
    
end

params{4} = amplitudes;
params{5} = intercepts;
end
