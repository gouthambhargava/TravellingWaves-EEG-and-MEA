function [allStableFrames] = getWaveStability(pgd,direction)
% Calculates the time points during which a wave can be said to be present.
% Takes into consideration the PGD values and the direction as well. 
% Output - is a vector of values that gives indices of the significant time
% points. Wave stability is calcualted from a zscored vector of values,
% where 0 has been set as the threshold. Can be varied to get either more
% time points with significant wave frames or less. 
% Inputs - the PGD and direction matrices (chan x time points) 

numElecs = numel(find(nansum(pgd,2))); %#ok<NANSUM>
stability = zeros(1,length(pgd));
for i = 1:size(pgd,1)
    for j = 1:size(pgd,2)-1
        waveStab = abs(pgd(i,j+1)*exp(sqrt(-1)*direction(i,j+1)) - pgd(i,j)*exp(sqrt(-1)*direction(i,j)));
        stability(j) = -sum(waveStab)/numElecs;
    end        
end
    allStableFrames = zscoreNan(stability);
    allStableFrames(allStableFrames<0) = 0;
end

function [nanScored] = zscoreNan(data) % zscoring when data has nans
    meanData = mean(data(~isnan(data)));
    stdData = std(data(~isnan(data)));
    nanScored = (data-meanData)/stdData;
end
    