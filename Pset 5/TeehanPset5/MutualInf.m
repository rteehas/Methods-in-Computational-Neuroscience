function mutual = MutualInf(probArray, P_T_n)

mutual = [];

s = size(probArray);

for i=1:s(1)
    theta_sum = 0;
    for j=1:s(2)
        n_sum = 0;
        for n=1:s(3)
            if probArray(i,j,n) ~= 0
                p = probArray(i,j,n) / P_T_n(i,n);
                pb = probArray(i,j,n)*log2(p);
                
            else
                pb = 0;
            end
            n_sum = n_sum + pb;
        end
%         p = probArray(i,j,:) ./ P_T_n(i,:);
%         n_sum = sum(probArray(i,j,:) .* log2(p), 2); 
%         theta_sum = theta_sum + ((1/13) * n_sum);
    theta_sum = theta_sum + ((1/13) * n_sum);
    end
    mutual = [mutual theta_sum];
end
            


end