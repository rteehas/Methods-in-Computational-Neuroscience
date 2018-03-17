function probArray = ConvertToArray(dirProbs, counts)

% get the max number of counts for the third dimension

z = max(max(max(counts)));

% indexed at 1 so add 1 to account for counts of 0

probArray = zeros(256,13,z+1);

for i=1:length(dirProbs)
    for j=1:length(dirProbs{1})
        cts = keys(dirProbs{i}{j});
        for k=1:length(cts)
            c = cts{k};
            probArray(j,i,c + 1) = dirProbs{i}{j}(c);
        end
    end
end


end