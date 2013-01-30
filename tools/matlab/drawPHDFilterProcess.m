close all;
clear all;

addpath '~/src/Matlab/export_fig' % export_fig can be downloaded from Matlab Central

drawLim = [-50, 50, -50, 50];
drawRes = 0.5;
x_min = drawLim(1);
x_max = drawLim(2);
y_min = drawLim(3);
y_max = drawLim(4);  

fig_width = 1600;
fig_height = 1200;

%% Birth Distribution

n = 5; % number of Gaussians

mean = zeros(n, 2);   % [x, y]
mean(1,:) = [20, 10];
mean(2,:) = [-30, 10];
mean(3,:) = [-15, -20];
mean(4,:) = [10, -10];
mean(5,:) = [-10, 30];

cov = zeros(n, 3);    % [xx, yy, xy]
cov(1,:) = [20, 20, 0];
cov(2,:) = [20, 20, 0];
cov(3,:) = [20, 20, 0];
cov(4,:) = [20, 20, 0];
cov(5,:) = [20, 20, 0];

weight = ones(n, 1); % [w]
weight(1) = 1;
weight(2) = 1;
weight(3) = 1;
weight(4) = 1;
weight(5) = 1;

figure;
[birth, h_birth] = drawGaussianMixture(mean, cov, weight, drawLim, drawRes);
ZLim([0, 0.02]);
colormap('summer')
set(gcf, 'Color', 'w');
set(gcf, 'Position', [0, 0, fig_width, fig_height]);
grid off;
set(gca, 'ZTick', []);
set(gca, 'ZTickLabel',[])
set(gca, 'box', 'off')
set(gca,'ZColor',get(gca,'Color'))

export_fig -r180 -opengl birth -pdf -eps -png

%% Prior Map

n = 3; % number of Gaussians

mean = zeros(n, 2);   % [x, y]
mean(1,:) = [25, -25];
mean(2,:) = [25, 35];
mean(3,:) = [-10, 25];

cov = zeros(n, 3);    % [xx, yy, xy]
cov(1,:) = [50, 50, 0];
cov(2,:) = [20, 20, 0];
cov(3,:) = [70, 20, -20];

weight = ones(n, 1); % [w]
weight(1) = 1;
weight(2) = 2;
weight(3) = 1.5;

figure;
[prior, h_prior] = drawGaussianMixture(mean, cov, weight, drawLim, drawRes);
ZLim([0, 0.02]);
colormap('summer')
set(gcf, 'Color', 'w');
set(gcf, 'Position', [0, 0, fig_width, fig_height]);
grid off;
set(gca, 'ZTick', []);
set(gca, 'ZTickLabel', [])
set(gca, 'box', 'off')
set(gca,'ZColor',get(gca,'Color'))

export_fig -r180 -opengl prior -pdf -eps -png

%% Prior Map merged with Birth

merged = birth + prior;
figure;
h_merged = surf(x_min : drawRes : x_max, y_min : drawRes : y_max, merged);
ZLim([0, 0.02]);
colormap('summer')
set(gcf, 'Color', 'w');
set(gcf, 'Position', [0, 0, fig_width, fig_height]);
grid off;
set(gca, 'ZTick', []);
set(gca, 'ZTickLabel', [])
set(gca, 'box', 'off')
set(gca,'ZColor',get(gca,'Color'))

export_fig -r180 -opengl merged -pdf -eps -png

%% Updated with Measurements

n = 15;

mean = zeros(n, 2);   % [x, y]
mean(1,:) = [20, 10];
mean(2,:) = [-30, 10];
mean(3,:) = [-15, -20];
mean(4,:) = [10, -10];
mean(5,:) = [-10, 30];
mean(6,:) = [25, -25];
mean(7,:) = [25, 35];
mean(8,:) = [-10, 25];
mean(9,:) = [-10, 0];
mean(10,:) = [-30, 30];
mean(11,:) = [40, -5];
mean(12,:) = [-30, 40];
mean(13,:) = [0, -35];
mean(14,:) = [-30, -10];
mean(15,:) = [40, -20];

cov = zeros(n, 3);    % [xx, yy, xy]
cov(1,:) = [20, 20, 0];
cov(2,:) = [20, 20, 0];
cov(3,:) = [20, 20, 0];
cov(4,:) = [20, 20, 0];
cov(5,:) = [20, 20, 0];
cov(6,:) = [50, 50, 0];
cov(7,:) = [20, 20, 0];
cov(8,:) = [70, 20, -20];
cov(9,:) = [30, 30, 0];
cov(10,:) = [40, 50, 10];
cov(11,:) = [30, 30, -10];
cov(12,:) = [60, 30, 0];
cov(13,:) = [30, 60, 20];
cov(14,:) = [45, 45, -20];
cov(15,:) = [25, 35, 0];

weight = ones(n, 1); % [w]
weight(1) = 1;
weight(2) = 1.2;
weight(3) = 0.5;
weight(4) = 1.5;
weight(5) = 1.6;
weight(6) = 1.7;
weight(7) = 2;
weight(8) = 1.5;
weight(9) = 0.1;
weight(10) = 0.2;
weight(11) = 0.3;
weight(12) = 0.4;
weight(13) = 0.5;
weight(14) = 0.6;
weight(15) = 0.4;

figure;
[update h_update] = drawGaussianMixture(mean, cov, weight, drawLim, drawRes);
ZLim([0, 0.02]);
colormap('summer')
set(gcf, 'Color', 'w');
set(gcf, 'Position', [0, 0, fig_width, fig_height]);
grid off;
set(gca, 'ZTick', []);
set(gca, 'ZTickLabel', [])
set(gca, 'box', 'off')
set(gca,'ZColor',get(gca,'Color'))

export_fig -r180 -opengl updated -pdf -eps -png

%% Pruned

pruneList = [9, 10, 11, 12, 13, 14, 15];
mean(pruneList,:) = [];
cov(pruneList,:) = [];
weight(pruneList,:) = [];

figure;
[pruned h_pruned] = drawGaussianMixture(mean, cov, weight, drawLim, drawRes);
ZLim([0, 0.02]);
colormap('summer')
set(gcf, 'Color', 'w');
set(gcf, 'Position', [0, 0, fig_width, fig_height]);
grid off;
set(gca, 'ZTick', []);
set(gca, 'ZTickLabel', [])
set(gca, 'box', 'off')
set(gca,'ZColor',get(gca,'Color'))

export_fig -r180 -opengl pruned -pdf -eps -png