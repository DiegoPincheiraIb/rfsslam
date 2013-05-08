function [likelihoodAllSeq] = multiFeatureLikelihood(T, t, c, Ec, showPermutations)

% Function multiFeatureLikelihood
% Multi-feature weighting approximation for the PHD Filter
% Keith Leung 2013
%
% Inputs:
%   T - likelihood table [number of measurements x number of features]
%   t - likelihood table threshold (entries below this value -> 0)
%   c - uniform clutter density
%   Ec - expected clutter
%
% Output:
%   likelihoodAllSeq

if nargin == 4
    showPermutations = 0;
end

likelihoodThreshold = t;
likelihoodAllSeq = 0;

% Feature sequence, interpret as z1 measures M(1), z2 measures M(2) ...
% If |Z| <= |M|, interpret A as Z, B as M
% If |Z| >  |M|, interpret A as M, B as Z
T_size = size(T);
if T_size(1) <= T_size(2);
    cardA = T_size(1);
    cardB = T_size(2);
else
    cardA = T_size(2);
    cardB = T_size(1);    
end
B = 1 : cardB;

isLastSequence = 0;
while(isLastSequence == 0)
    
    % Find the likelihood of the current sequence
    likelihoodCurSeq = 1;

    for i = 1:cardA
        
        if(T_size(1) <= T_size(2))
            likelihoodCurPair = T(i, B(i));
        else
            likelihoodCurPair = T(B(i), i);
        end
            
        % Is likelihood pair for measurement i and feature M(i) small?
        if likelihoodThreshold <= 0 || likelihoodCurPair > likelihoodThreshold
            % Measurement likelihood is not small
            likelihoodCurSeq = likelihoodCurSeq * likelihoodCurPair;

        else
            % Measurement likelihood is small
            likelihoodCurSeq = 0;

            % There is no point consider sequences with M(i)
            % Fast foward to the last sequence that has M(i)
            if i < cardB
                B(i+1:cardB) = sort(B(i+1:cardB), 'descend');
            end
            break;

        end

    end
    
    %if(showPermutations == 1 && likelihoodCurSeq ~= 0)
    %    disp(B);
    %end
    %disp(B);

    likelihoodAllSeq = likelihoodAllSeq + likelihoodCurSeq;
    
    % If there are more features than measurements
    % M(cardZ + 1 : cardM) are not used, so fast forward
    if cardA < cardB
        B(cardA + 1 : cardB) = sort( B(cardA + 1 : cardB), 'descend');
    end
    
    % Find the next sequence in lexicographical order
    for i = cardB-1 : -1 : 0

        if i == 0
            isLastSequence = 1;
            break
        end

        if B(i) < B(i+1)

            % Swap M(i) with the highest index j, such that M(j) > M(i)
            for j = cardB: -1 : 0

                if B(j) > B(i)

                    temp = B(j);
                    B(j) = B(i);
                    B(i) = temp;
                    break;

                end

            end

            % Reverse order of elements after M(i), denote as M_
            cardB_ = cardB - i;
            for j = 1 : floor(cardB_/2)

                swapIdx1 = i + j;
                swapIdx2 = cardB - j + 1;
                temp = B( swapIdx2 );
                B(swapIdx2) = B(swapIdx1);
                B(swapIdx1) = temp;

            end
            
            break;

        end

    end

    % At this point we have got the next sequence

end

% Account for clutter if |Z| > |M|
if T_size(1) > T_size(2)
    likelihoodAllSeq = likelihoodAllSeq * c^(T_size(1) - T_size(2)) / Ec;
end

