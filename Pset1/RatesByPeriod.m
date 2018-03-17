function rates = RatesByPeriod(trials)

pretrial = [];
sample = [];
delay = [];
test = [];

for i=1:length(trials)
    t = trials{i};
    pre = t(t < 0);
    samp = t(t >= 0 & t <= .65);
    del = t(t >.65 & t <= 1.65);
    tes = t(t > 1.65);
    
    if isempty(pre)
        pre_rate = 0;
    else
        pre_rate = length(pre) / (-pre(1));
    end
    if isempty(samp)
        samp_rate = 0;
    else
        samp_rate = length(samp) / (samp(length(samp)) - samp(1));
    end
    
    if isempty(del)
        del_rate = 0;
    else
        del_rate = length(del) / (del(length(del)) - del(1));
    end
    
    if isempty(tes)
        tes_rate = 0
    else
        tes_rate = length(tes) / tes(length(tes));
    end
    
    pretrial = [pretrial pre_rate];
    sample = [sample samp_rate];
    delay = [delay del_rate];
    test = [test tes_rate];
  
end

rates = {pretrial, sample, delay, test};
end