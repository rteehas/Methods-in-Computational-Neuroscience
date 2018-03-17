function [y,t] = Current(t, A)

% frequency of 40hz
freq = 40;

% we want the function to be entirely positive, so we add A/2 
% and have the amplitude of the sine function be A/2

y = A*sin(2 * pi * freq * t);


% only select when current is positive or 0
t = t(y>=0);

y = y(y>=0);
end