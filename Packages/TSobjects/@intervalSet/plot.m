function plot(is)
% makes a plot to represent an intervalset but simply plotting horizontal
% lines based on the start and stop time of each interval.  Units of plot
% is in seconds.
%
% INPUTS:
% is - an intervalSet
%
% OUTPUTS:
% (only graphical)

% Brendon Watson 2014

SampFreq = 10000;

spans = [Start(is) End(is)];
spans = spans/SampFreq;

plot(spans',zeros(size(spans')))