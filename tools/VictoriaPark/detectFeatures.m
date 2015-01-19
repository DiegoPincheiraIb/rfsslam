% Use detectTreesI16.m to detect features in aa3.lsr2.mat for all timesteps

format long
clear all
load aa3_lsr2.mat

% Fix timestamps

TLsr = double(TLsr - TLsr(1) + 852) / 1000;

% Extract features 
global AAr; 
AAr = [0:360]*pi/360 ;
detections = zeros( length(TLsr) * 20, 4 );
nDetections = 0;
mask13 = uint16(2^13 -1) ;

for t = 1 : length(LASER(:,1))
   
    scan = double( bitand( mask13,LASER(t,:) ) ) / 100;
    trees = detectTrees(scan);
    if isempty(trees)
        continue;
    end
    nTrees = length(trees(1,:));
    for m = 1 : nTrees
        nDetections = nDetections + 1;
        detections( nDetections, 1 ) = TLsr(t);
        detections( nDetections, 2:4 ) = trees(:,m)';
    end
    
end

detections(nDetections + 1 : end, :) = [];

disp('Max range')
disp( max(detections(:,2)) )
disp('Max diameter')
disp( max(detections(:,4)) )
disp('Min diameter')
disp( min(detections(:,4)) )

% Write to file
FID = fopen('measurements.dat','w');
fprintf(FID, '%10.3f %10.5f %10.5f %10.5f\n', detections');
fclose(FID);