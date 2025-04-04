% sPos: position of sf or size parameter. Takes values between 1-6: [0.5 1 2 4 8 all]
% oriPos: position of orientation parameter. Takes values between 1-9: [0 22.5 45 67.5 90 112.5 135 157.5 all] 
% param: 'sf' or 'size' to specify whether this corresponds to an SFOri or a sizeOri protocol

function [allData,goodElectrodes,timeVals,rfData,parameters] = loadData(subjectName,expDate,protocolName,dataPath,gridType,sPos,oriPos,param)

if ~exist('param','var');       param = 'sf';                           end

folderName = fullfile(dataPath,subjectName,gridType,expDate,protocolName);

% Get good electrodes    
rfData = load([subjectName gridType 'RFData.mat']);
goodElectrodes = rfData.highRMSElectrodes;
goodElectrodes = goodElectrodes(goodElectrodes<=81); % Only microelectrodes
numGoodElectrodes = length(goodElectrodes);

% Get good trials
parameters = load(fullfile(folderName,'extractedData','parameterCombinations.mat'));
t = load(fullfile(folderName,'segmentedData','LFP','lfpInfo.mat'));
timeVals = t.timeVals;


badTrials = load(fullfile(folderName,'segmentedData','badTrials.mat'));
badTrials = badTrials.badTrials;
if strcmp(param,'sf')
    goodPos = setdiff(parameters.parameterCombinations{1,1,1,sPos,oriPos},badTrials);
elseif strcmp(param,'size')
    goodPos = setdiff(parameters.parameterCombinations{1,1,sPos,1,oriPos},badTrials);
end

allData = zeros(numGoodElectrodes,length(goodPos),length(timeVals));

for i=1:numGoodElectrodes
    lfpData = load(fullfile(folderName,'segmentedData','LFP',['elec' num2str(goodElectrodes(i)) '.mat']));
    allData(i,:,:) = lfpData.analogData(goodPos,:);
end
end