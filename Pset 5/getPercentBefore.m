function frac = getPercentBefore(m, times, cutoff)


% record cutoff after neural response, so after it is most rapidly
% increasing
[val argmax] = max(diff(m));

pt = times(argmax) + cutoff;

frac = m(times == pt) / m(256);


end