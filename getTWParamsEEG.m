function outputs = getTWParamsEEG(phiMat,burstMat,timeVals,locList,boundryLims,badElecs)
% Inputs
% phiMat - phases of single trial data in electrodes x time points format
% burstMat - time series indicting burst times in electrodes x time points format
% timeVals - vector of time values
% goodElectrodes - list of good electrodes
% locList - array of size electrodes x 2 containing the x and y position of each electrode
% boundryLims - defines the timepoints of interest while identifying waves.Ex. [0.25 0.75]
% badElecs - list of bad electrodes that can be supressed during wave
% calculation
% nPerm - when not set to [], TW statistics are computed nPerm times after shuffling the phases.

% if ~exist('nPerm','var'); nPerm = [];  end
if ~exist('badElecs','var'); badElecs = []; end 

% Parameters
% fs = 1/(timeVals(2)-timeVals(1)); %sampling frequency
timePoints = length(timeVals);
goodElectrodes = 1:size(phiMat,1);
numGoodElectrodes = length(goodElectrodes);
boundryLims = dsearchn(timeVals',boundryLims(1)):dsearchn(timeVals',boundryLims(2));

%initialize results
pgd = nan(numGoodElectrodes,timePoints);
direction = nan(numGoodElectrodes,timePoints);
sFreq = nan(numGoodElectrodes,timePoints);

burstMat(isnan(burstMat)) = 0;
burstMat(badElecs,:) = 0;

%% run circ reg after getting significant electrodes from burst detection
tic
for timei = 1:numel(boundryLims)
    phaseData = phiMat(:,boundryLims(timei));
    burstElecs = burstMat(:,boundryLims(timei));
    % do regression analysis on the polar and linear coordinates
    [directionT, pgdT,spatialFreqT] = getWaveClusterAnalysis(locList,phaseData,burstElecs,2);
    direction(:,boundryLims(timei)) = directionT;
    pgd(:,boundryLims(timei)) = pgdT;
    sFreq(:,boundryLims(timei)) = spatialFreqT;
end
toc
% clean up the pgd and direction outputs
pgd(pgd<=0) = nan; % remove nagative PGD values
direction(isnan(pgd)) = nan; % so that directions show up only during significant PGD points.

%% find frames of tiem series that showed a good wave
allStableFrames = getWaveStability(pgd,direction);

%%  get additional params
outputs.direction = direction+pi;
outputs.pgd = pgd;
outputs.waveFrames = allStableFrames;
outputs.sFreq = sFreq;
deltaT = 1/2000; % in s
phi_dot = circ_mean(circ_dist(phiMat(:,2:end),phiMat(:,1:end-1)));
tempFreq = (phi_dot/deltaT)/(2*pi); % convert to Hz
outputs.speed = ([tempFreq,0]./sFreq)/1000; % convert to m/s
outputs.speed(outputs.speed==inf) = nan;
end