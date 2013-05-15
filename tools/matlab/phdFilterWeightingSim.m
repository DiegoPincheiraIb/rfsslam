% PHD Filter Weighting Sim
% Keith Leung 2013

% This script is for testing different approximations for the weighting
% strategy in the PHD Filter, with focus on the measurement likelihood term

clear all
close all

% Settings
workspace = [-10, -10, 10, 10];
robot = [0; 0];
nFeatures_max = 7;
nParticles = 1;
R = [0.1,0; 0,0.1];
P_detection = 0.95;
applyThreshold = 0;
likelihoodThreshold = 0.05;
clutterExpected = [1.0 5.0 20.0];
clutterExpectedExp = exp(clutterExpected);
clutterDensity = clutterExpected / 20^2;

Pd_trial = 1:-0.025:0.1;

ml_sums = zeros(nFeatures_max, length(Pd_trial), length(clutterExpected));
ml_errors = zeros(nFeatures_max, length(Pd_trial), length(clutterExpected));

for nFeatures = 1:nFeatures_max
for Pd_trial_counter = 1:length(Pd_trial)
for clutterCounter = 1:length(clutterExpected);
    
    clutterExpectedExp = exp(clutterExpected(clutterCounter));
    clutterDensity = clutterExpected(clutterCounter) / 20^2;
    
    
P_detection = Pd_trial(Pd_trial_counter);
display(sprintf('trial %d -- %f -- %f', nFeatures, P_detection, clutterExpected(clutterCounter)));

% Randomly generate a bunch of features in 2d, also to serve as prior map
map = [workspace(3) - workspace(1), 0; 0, workspace(4) - workspace(2)]*rand(2, nFeatures) +  [workspace(1); workspace(2)]*ones(1, nFeatures);
S = [0.5, 0; 0, 0.5];

% Generate measurements
z = map - robot * ones(1, nFeatures);
nObs = nFeatures;

% Generate particles
p = [workspace(3) - workspace(1), 0; 0, workspace(4) - workspace(2)]*rand(2, nParticles) +  [workspace(1); workspace(2)]*ones(1, nParticles);
n = ceil(nParticles / 5); % 20% of particles close to true position
p(:, 1:n) = robot * ones(1,n) + randn(2, n);
p(:, 1) = robot; % first particle is at true position

% Weighting
pTable = perms(1:nFeatures);
Gaussian_multiplier = ((2*pi)^2 * det(R))^0.5;
for i = 1:nParticles
    mxzLikelihood = zeros(nFeatures, nObs);
    
    % Expected measurements
    z_exp = map - p(:,i)*ones(1, nFeatures);    
    
    for m = 1:nFeatures
        
        for nz = 1:nObs;
            
            % Measurement likelihood table
            mDist = (z(:,nz) - z_exp(:,m))' / R * (z(:,nz) - z_exp(:,m));
            mxzLikelihood(m, nz) = Gaussian_multiplier * exp(-0.5 * mDist);
            
            if applyThreshold == 1 && mxzLikelihood(m, nz) < likelihoodThreshold
                mxzLikelihood(m, nz) = 0;
            end
        
        end
    
    end
    
    % Set measurement likelihood
    ml = zeros(nObs + 1, 1);
    for n_features_observed = 0:nObs
        
        n_features_not_observed = nObs - n_features_observed;
        
        f1 = (1 - P_detection)^(n_features_not_observed); % Prob of missed detection
        f2 = (P_detection)^(n_features_observed); % Prob of detection
        f3 = clutterDensity^(nObs - n_features_observed) / clutterExpectedExp; % clutter
        f4 = 0; % sum of likelihoods (to be updated)
        
        if applyThreshold == 0
            if n_features_observed == 0
                f4 = 1;
            else 
                featureCombinations = nchoosek(1:nFeatures, n_features_observed);
                measurementCombinations = nchoosek(1:nObs, n_features_observed);
                measurementPermutationTable = pTable(1:factorial(n_features_observed), n_features_not_observed + 1 : end); 

                % Pick a combination of feature indices
                for featureCombo = 1 : length( featureCombinations(:,1) )

                    featureCombination = featureCombinations( featureCombo, : );

                    % Pick a combination of measurements
                    for obsCombo = 1 : length( measurementCombinations(:,1) )

                        observationCombination = measurementCombinations( obsCombo, : );

                        % Pick a permutation of measurements
                        for obsPermu = 1 : length( measurementPermutationTable(:,1) )

                            observationPermutation = observationCombination( measurementPermutationTable( obsPermu, : ) );
                            %display( [featureCombination, observationPermutation] );

                            likelihoodProduct = 1;
                            for i = 1:length( observationPermutation )

                                m_idx = featureCombination(i);
                                z_idx = observationPermutation(i);
                                likelihoodProduct = likelihoodProduct * mxzLikelihood(m_idx, z_idx);

                            end
                            f4 = f4 + likelihoodProduct;

                        end

                    end

                end

            end
            
        elseif applyThreshold == 1
            
            f4 = 0;
            likelihoodProducts = zeros(factorial(1, nFeatures));
            likelihoodProducts(1:length( mxzLikelihood(1,:) )) = mxzLikelihood(1,:);
            
            for i = 2:length(mxzLikelihood(:,1)) % row
                
                
            end
            
            f4 = f4 + likelihoodProdclutterExpected
            
            
        end
        
        ml(n_features_observed + 1) = f1 * f2 * f3 * f4;
        
    end
    
end

ml_sums(nFeatures, Pd_trial_counter, clutterCounter) = sum(ml(1:n_features_observed + 1));
ml_errors(nFeatures, Pd_trial_counter, clutterCounter) = sum(ml(1:n_features_observed));
%display(ml_sum);

end
end
end

display(ml_errors)
save('approximationError_0.2Clutter.mat');

figure;
plot_nFeature = 5;
for c = 1:length(clutterExpected)
    plot( Pd_trial, ml_errors(plot_nFeature, :, c) ./ ml_sums(plot_nFeature, :, c) * 100, '-k')
    hold on
end
%set(gca,'YScale','log');
%title('Measurement Likelihood Error');
xlabel('Probability of detection');
ylabel('Percantage Error');
grid on
set(gcf, 'Color', 'w');

addpath '~/src/Matlab/export_fig' % export_fig can be downloaded from Matlab Central
export_fig(sprintf('results/likelihoodError_%f_thresholded.pdf',clutterExpected));

% % Felipe's method
% 
% f = P_detection ^ nFeatures;
% zSums = sum(mxzLikelihood);
% zProd = 1;
% for i = 1:length(zSums)
%     zProd = zProd * zSums(1); 
% end
% display(f * zProd /  nClutterExpected);
% 
% 
% plot(map(1,:), map(2,:), 'r.');
% xlim([workspace(1), workspace(3)]);
% ylim([workspace(2), workspace(4)]);
% axis square;
% hold on
% plot(p(1,:), p(2,:), 'b.');