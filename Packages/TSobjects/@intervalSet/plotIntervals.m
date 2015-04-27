function plotIntervals(intervals)
% % Plots a horizontal line showing timing of each interval in a set
%
% USAGE:
%     plotIntervals(x)
%
% INPUT:
%    x: intervalSet array to be plotted
%
% Brendon Watson, 2014


sta = Start(intervals);
sto = End(intervals);

plot([sta sto]',zeros(size([sta sto]')))