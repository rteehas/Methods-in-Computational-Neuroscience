function P_T_n = getCountProbs(probArray)

P_T_n = zeros(256, 24);

for i=1:length(probArray)
    for j=1:length(probArray(1,1,:))
        P_T_n(i,j) = sum(probArray(i,:,j) * (1/13), 2);
    end
    
end

end