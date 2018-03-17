function distances = GetDistances(i,j,k,spikes, q)
distances = {};
distances{i} = {};
distances{i}{j} = zeros(length(spikes{i}), length(spikes{i}{j}));
train = spikes{i}{j}{k};

for i1=1:length(spikes{1})
    for j1=1:length(spikes{1}{1})
           distances{i}{j}(i1,j1) = spkd(train, spikes{i}{i1}{j1}, q);
            
    end
end


end