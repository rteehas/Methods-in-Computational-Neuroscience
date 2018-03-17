function probs = GetProbs(data, dir_index, directions)



probs = {};

d = directions(dir_index);

for i=1:length(data)
    probs{i} = containers.Map('KeyType','double','ValueType','double');
    for j=1:length(data(i, dir_index, :))
        n = data(i, dir_index, j);
        if ~probs{i}.isKey(n)
            probs{i}(n) = 1;
        else 
            probs{i}(n) = probs{i}(n) + 1;
        end
        
    end
    
end


for i=1:length(probs)
    k = keys(probs{i});
    for j=1:length(k)
        ct = k{j};
        probs{i}(ct) = probs{i}(ct) / 184;
    end
end






end