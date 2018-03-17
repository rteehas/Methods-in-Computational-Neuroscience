function [x_new, y_new] = selectTime(x,y,spikes)

x_new = x(spikes ==1);
y_new = y(spikes ==1);


end