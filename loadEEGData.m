function [eegData,badElectrodes,timeVals] = loadEEGData(subjectIndex,dataPath)

% load data for TLSA
projectName = 'ADGammaProject';
% protocolType = 'SF_ORI';
% folderName = 'G:\TLSAdata';
dataType = 'cleanData';

% get list of good subjects
goodSubjectsList = getGoodSubjectsProjectwise(projectName,1,protocolType);
goodSubjectsList = goodSubjectsList{1};

% get exp details of all good subjects
for i = 1:numel(goodSubjectsList)
    [expDates,protocolNames,capLayout,usableDataFlag(i)] = getProtocolDetailsForAnalysis(projectName,goodSubjectsList{i},protocolType);
    if ~isempty(expDates)
        expDatesAll{i} = expDates{1};
        protocolNamesAll{i} = protocolNames{1};
        capLayoutAll{i} = capLayout{1};
    end
end

if ~isempty(protocolNamesAll{subjectIndex})
    subData = load(fullfile(dataPath,dataType,protocolType,strcat(goodSubjectsList{subjectIndex},'-',expDatesAll{subjectIndex},'-',protocolNamesAll{subjectIndex},'.mat')));
    eegData = subData.eegData;
    timeVals = subData.timeVals;
    badElectrodes = subData.badElecs.badImpedanceElecs;
    % goodElectrodes = setdiff(1:size(eegData,1),goodElectrodes);
end

end












