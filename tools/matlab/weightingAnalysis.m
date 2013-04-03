% RB-PHD-FILTER weighting analysis
% Keith Leung 2013
clear all
close all

for c = 0:0.1:2
    x = 1:100;
    for i = 1:100
        y(i) = (x(i) / (c + x(i)));
    end
    plot(x,y)
    hold on
end

grid on