function [direction, pgd,spatialFreq] = getWaveClusterAnalysis(elecLocs,phaseData,burstLocs,testCase)

% Inputs
% elecLocs is the list of electrode locations (x and y locations in the grid - dimensions reduced from 3D)
% phaseData is the instant phases of all elecs at the time point
% burstLocs is a logical vector, same length as the phaseData with gamma burst locations 
% testCase - Specify method being used to identify electrodes to be
% processed
%     1 - uses all electrodes at the time point that have a burst. Generates a single value for all electrodes.
%     2 - uses the distance matrix generated electrode locations. Only calculates direction and PGD for the selected electrode and the nearest 3 neighbours. 
%         Provides a unique direction value for each electode. No of neighbours can be adjusted. 
% Outputs
% Also described in circRegressEEG. 
%%        
% initialize the outputs
direction = nan(1,length(phaseData));
pgd = nan(1,length(phaseData));
spatialFreq = nan(1,length(phaseData));

% distmat = load('D:\IISC_work\gitScripts\wavePatternAnalysis\distMat.mat');
% distmat = distmat.distMat;
distmat = squareform(pdist(elecLocs)); % define the distance matrix
burstElecs = find(burstLocs); %find indices of electrodes that have bursts
numElecs = 1:length(phaseData); % find total number of electrodes

if ~isempty(burstElecs)
    if testCase == 1
        % method 1 - consider all electrodes
        newPhases = phaseData(burstElecs);
        newLocs = elecLocs(burstElecs,:);        
        [directionI,spatialFreqI,pgdI] = circRegressEEG(newPhases,newLocs);
        direction(burstElecs) = directionI;
        spatialFreq(burstElecs) = spatialFreqI;
        pgd(burstElecs) = pgdI;    
    elseif testCase == 2 
        % method 2 - select electrode clusters based on distance
        for i = 1:length(elecLocs)
            [~,sortedIndices] = sort(distmat(i,:)); % sort the electrodes based on distance from current electrode in the loop
            % neighbourElecs = sortedIndices(1:4);
            sortedIndices(setdiff(numElecs,burstElecs)) = [];
            sigElecs = intersect(neighbourElecs,burstElecs);
            newLocs = elecLocs(sigElecs,:);
            newPhases = phaseData(sigElecs,:);
            if numel(newPhases)>2 
                [direction(i),spatialFreq(i),pgd(i)] = circRegressEEG(newPhases,newLocs);
            end
        end
    end
end