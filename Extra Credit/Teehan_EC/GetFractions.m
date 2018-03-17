function frac = GetFractions(spikes, qs)


frac = [];
    
for l=1:length(qs)
    total = 0;
    true = 0;
    for i=1:length(spikes)
        for j0=1:length(spikes{1})
            for k=1:length(spikes{1}{1})
                    d = GetDistances(i,j0,k, spikes,qs(l));
                    [val arg] = min(mean(d{i}{j0},2));
                    if arg == j0
                        true = true + 1;
                    end
                    total = total+1;

            end
        end
    end
    frac = [frac true/total];
end


end