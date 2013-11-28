clear all
close all

%logDir = '../../data/rbphdslam/';
%logDir = '../../data/fastslam/';
logDir = '../../data/mhfastslam/';

fid_mapError = fopen(strcat(logDir, 'landmarkEstError.dat'));

mapError = textscan(fid_mapError, '%d %d %f %f\n', [4 inf]);
timesteps = mapError{1}(:);
ospa_error = mapError{4}(:);
map_count = mapError{2}(:);
est_count = mapError{3}(:);

figure;
plot(timesteps, ospa_error,'r-');
ylim([0 3]);
xlabel('Timestep');
ylabel('OSPA Error');
%export_fig /home/kykleung/Desktop/rbphdslam_M8_pd99c01_OSPAError.pdf
%export_fig /home/kykleung/Desktop/fastslam_H1_pd99c01_OSPAError.pdf
%export_fig /home/kykleung/Desktop/mhfastslam_H1_pd99c01_OSPAError.pdf
export_fig(strcat(logDir, 'OSPAError.pdf'));

figure;
plot(timesteps, map_count,'k-','LineWidth',2);
hold on
plot(timesteps, est_count,'r-');
xlabel('Timestep');
ylabel('Landmark Count');
ylim([0 70]);
legend('Location', 'SouthEast', 'Actual', 'Estimated');
%export_fig /home/kykleung/Desktop/rbphdslam_M8_pd99c01_cardinalityError.pdf
%export_fig /home/kykleung/Desktop/fastslam_H1_pd99c01_cardinalityError.pdf
%export_fig /home/kykleung/Desktop/mhfastslam_H1_pd99c01_cardinalityError.pdf
export_fig(strcat(logDir, 'cardinalityError.pdf'));