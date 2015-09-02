function h = plot(S)
% Plots a tsdArray as raster (based on time only, no values used).
% Brendon Watson 2015

x = [];
y = [];
S = cellArray(S);

for a = 1:length(S);
    tx = Range(S{a},'s');
    ty  = a*ones(size(tx));
    x = cat(1,x,tx);
    y = cat(1,y,ty);
end

plot(x,y,'.')
axis tight