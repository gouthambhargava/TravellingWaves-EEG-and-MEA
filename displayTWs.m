function displayTWs(subjectName,expDate,protocolName,dataPath,sPos,oriPos,selectedElectrodes,freqRangeList,dataType)

% Input parameters
% subjectName - for monkey data define subject name. Otherwise define subject index for EEG data
% expData - for monkey data only
% protocolName - define protocolName - only for monkey data
% sPos and oriPos - the index for spatial frequency and orientation - for monkey data only
% selectedElectrodes - define 3 electrode indices to display
% freqRangeList - define the required frequency ranges. No of plots depends on number of ranges defined. 
% dataType - 1 for monkey data (MEA) and 2 for EEG data.

fontSizeSmall = 10; fontSizeMedium = 12; fontSizeLarge = 16;
backgroundColor = 'w'; panelHeight = 0.125;

numFrequencyRanges = length(freqRangeList);
colorNamesFreqRanges = gray(numFrequencyRanges+1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Load data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Getting data...');
if dataType ==1
    gridType = 'Microelectrode';  
    [allData,goodElectrodes,timeVals,rfData,parameters] = loadMData(subjectName,expDate,protocolName,dataPath,gridType,sPos,oriPos);
    numGoodElectrodes = length(goodElectrodes);
    % get location list
    locList = nan(numGoodElectrodes,2);
    for iElec = 1:numGoodElectrodes
        [locList(iElec,1),locList(iElec,2)] = find(electrodeArray==goodElectrodes(iElec));
    end
else 
    [allData,badElectrodes,timeVals] = loadEEGData(subjectIndex,dataPath);
    goodElectrodes = setdiff(1:size(allData,1),badElectrodes);
    numElectrodes = size(allData,1);
    % generate the location list
    load('D:\IISC_work\gitScripts\Dependancies\Programs\ProgramsMAP\Montages\actiCap64.mat')
    locList = zeros(numElectrodes,2);
    elecInfo = [chanlocs(1:numElectrodes).theta;chanlocs(1:numElectrodes).radius];
    theta = elecInfo(1,:)';
    radius = elecInfo(2,:)';
    [locList(:,1), locList(:,2)] = pol2cart((pi/180).*((-1).*theta+90),radius); % X, Y channel coords
end
disp('Data loaded')
%%%%%%%%%%%%%%%%%%%%%%% Do burst estimation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% shift this part further down and add more code so that other methods of
% burst detection are also determined - including modal and none.
% thresholdFactor = 4;
% stimulusDurationS = [0 0.8]; % Stimulus duration to be highlighted
% baselinePeriodS = [-0.5 0];
analysisWindow = [0.25 0.75];
% analysisPeriodS = [-0.5 1];
% filterOrder = 4;

disp('Performing burst analysis...');
burstTS = zeros(size(allData,1),size(allData,3),numFrequencyRanges);
filteredSignal = zeros(size(allData,1),size(allData,3),numFrequencyRanges);
bandPhase = zeros(size(allData,1),size(allData,3),numFrequencyRanges);
bandPower = zeros(size(allData,1),size(allData,3),numFrequencyRanges);
for iFreq = 1:numFrequencyRanges
    if strcmp(analysisMethod,'hilbert')
        [burstTS(:,:,iFreq),filteredSignal(:,:,iFreq),bandPhase(:,:,iFreq), bandPower(:,:,iFreq)] = filterTW(allData,srate,freqRangeList{iFreq},analysisWindow,2,timeVals);
    elseif strcmp(analysisMethod,'MODAL')
        [burstTS(:,:,iFreq),filteredSignal(:,:,iFreq),bandPhase(:,:,iFreq), bandPower(:,:,iFreq)] = filterTW(allData,srate,freqRangeList{iFreq},analysisWindow,1,timeVals);
    else
        [burstTS(:,:,iFreq),filteredSignal(:,:,iFreq),bandPhase(:,:,iFreq), bandPower(:,:,iFreq)] = filterTW(allData,srate,freqRangeList{iFreq},analysisWindow,3,timeVals);
    end
end
disp('Burst detection complete')

%%%%%%%%%%%%%%%%%%%%%%%%%%%  Data selection %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hPanel1 = uipanel('Title','Data','fontSize',fontSizeLarge,'Unit','Normalized','Position',[0 1-panelHeight 0.2 panelHeight]);

% Selected electrodes
uicontrol('Parent',hPanel1,'Unit','Normalized','Position',[0 0.5 0.25 0.5],'Style','text','String','Elecs','FontSize',fontSizeMedium);
hElectrodes = uicontrol('Parent',hPanel1,'Unit','Normalized','BackgroundColor', backgroundColor, 'Position', [0.25 0.5 0.75 0.5],'Style','edit','String',num2str(selectedElectrodes),'FontSize',fontSizeMedium);

% Trial Number
numTrials = size(allData,2);
uicontrol('Parent',hPanel1,'Unit','Normalized','Position',[0 0 0.25 0.5],'Style','text','String','TrialNum','FontSize',fontSizeMedium);
trialNumList = 1:numTrials;
hTrial = uicontrol('Parent',hPanel1,'Unit','Normalized','BackgroundColor', backgroundColor, 'Position', [0.25 0 0.15 0.5],'Style','popup','String',trialNumList,'FontSize',fontSizeMedium);

% Plot data
uicontrol('Parent',hPanel1,'Unit','Normalized','Position',[0.4 0 0.2 0.5],'Style','pushbutton','String','plot','FontSize',fontSizeMedium,'Callback',{@plot_Callback});
uicontrol('Parent',hPanel1,'Unit','Normalized','Position',[0.6 0 0.2 0.5],'Style','pushbutton','String','Rescale','FontSize',fontSizeMedium,'Callback',{@rescale_Callback});
uicontrol('Parent',hPanel1,'Unit','Normalized','Position',[0.8 0 0.2 0.5],'Style','pushbutton','String','Clear','FontSize',fontSizeMedium,'Callback',{@cla_Callback});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Axis Ranges %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hPanel2 = uipanel('Title','AxisRanges1','fontSize',fontSizeLarge,'Unit','Normalized','Position',[0.2 1-panelHeight 0.2 panelHeight]);
axisRange1List0{1} = [0 80]; axisRange1Name{1} = 'FreqLims (Hz)';
axisRange1List0{2} = [-0.5 1]; axisRange1Name{2} = 'TimeLims (S)';
axisRange1List0{3} = [-2 2]; axisRange1Name{3} = 'cLims (raw)';

numAxisRanges1 = length(axisRange1List0);
hAxisRange1Min = cell(1,numAxisRanges1);
hAxisRange1Max = cell(1,numAxisRanges1);

for ii=1:numAxisRanges1
    uicontrol('Parent',hPanel2,'Unit','Normalized','Position',[0 1-ii/numAxisRanges1 0.5 1/numAxisRanges1],'Style','text','String',axisRange1Name{ii},'FontSize',fontSizeSmall);
    hAxisRange1Min{ii} = uicontrol('Parent',hPanel2,'Unit','Normalized','BackgroundColor', backgroundColor,'Position',[0.5 1-ii/numAxisRanges1 0.25 1/numAxisRanges1], ...
        'Style','edit','String',num2str(axisRange1List0{ii}(1)),'FontSize',fontSizeSmall);
    hAxisRange1Max{ii} = uicontrol('Parent',hPanel2,'Unit','Normalized','BackgroundColor', backgroundColor,'Position',[0.75 1-ii/numAxisRanges1 0.25 1/numAxisRanges1], ...
        'Style','edit','String',num2str(axisRange1List0{ii}(2)),'FontSize',fontSizeSmall);
end

hPanel3 = uipanel('Title','AxisRanges2','fontSize',fontSizeLarge,'Unit','Normalized','Position',[0.4 1-panelHeight 0.2 panelHeight]);

axisRange2List0 = cell(1,numFrequencyRanges);
axisRange2Name = cell(1,numFrequencyRanges);
for iFreq=1:numFrequencyRanges
    axisRange2List0{iFreq} = [-50 50];
    axisRange2Name{iFreq} = ['yRange (' num2str(freqRangeList{iFreq}(1)) '-' num2str(freqRangeList{iFreq}(2)) ' Hz)'];
end

hAxisRange2Min = cell(1,numFrequencyRanges);
hAxisRange2Max = cell(1,numFrequencyRanges);

for ii=1:numFrequencyRanges
    uicontrol('Parent',hPanel3,'Unit','Normalized','Position',[0 1-ii/numFrequencyRanges 0.5 1/numFrequencyRanges],'Style','text','String',axisRange2Name{ii},'FontSize',fontSizeSmall);
    hAxisRange2Min{ii} = uicontrol('Parent',hPanel3,'Unit','Normalized','BackgroundColor', backgroundColor,'Position',[0.5 1-ii/numFrequencyRanges 0.25 1/numFrequencyRanges], ...
        'Style','edit','String',num2str(axisRange2List0{ii}(1)),'FontSize',fontSizeSmall);
    hAxisRange2Max{ii} = uicontrol('Parent',hPanel3,'Unit','Normalized','BackgroundColor', backgroundColor,'Position',[0.75 1-ii/numFrequencyRanges 0.25 1/numFrequencyRanges], ...
        'Style','edit','String',num2str(axisRange2List0{ii}(2)),'FontSize',fontSizeSmall);
end

%%%%%%%%%%%%%%%%%%%%  Traveling Wave plot panel %%%%%%%%%%%%%%%%%%%%%%%%
%this needs to be still be implemented
hPanel4 = uipanel('Title','Traveling Wave','fontSize',fontSizeLarge,'Unit','Normalized','Position',[0.6 1-panelHeight 0.2 panelHeight]);
uicontrol('Parent',hPanel4,'Unit','Normalized','Position',[0 0.5 0.25 0.5],'Style','text','String','Electrode Fraction','FontSize',fontSizeSmall);
hElecFrac = uicontrol('Parent',hPanel4,'Unit','Normalized','BackgroundColor', backgroundColor,'Position',[0.25 0.65 0.2 0.3], ...
    'Style','edit','String','0.6','FontSize',fontSizeSmall);

uicontrol('Parent',hPanel4,'Unit','Normalized','Position',[0 0 0.25 0.5],'Style','text','String','Select Electrodes','FontSize',fontSizeSmall);
hSelectElec = uicontrol('Parent',hPanel4,'Unit','Normalized','BackgroundColor', backgroundColor,'Position',[0.25 0.2 0.2 0.3], ...
    'Style','popup','String',{'All','Only Immediate neighbours','Extended neighbours'},'FontSize',fontSizeSmall);

uicontrol('Parent',hPanel4,'Unit','Normalized','Position',[0.5 0.55 0.25 0.45],'Style','text','String','Ref Choice','FontSize',fontSizeSmall);
hSelectBurstMet = uicontrol('Parent',hPanel4,'Unit','Normalized','BackgroundColor', backgroundColor,'Position',[0.75 0.72 0.2 0.3], ...
    'Style','popup','String',{'Hilbert','MODAL','None'},'FontSize',fontSizeSmall);

% these need a callback function
plotToggleTW1 = uicontrol('Parent',hPanel4,'Unit','Normalized','Position',[0.5 0.2 0.25 0.3],'Style','togglebutton','String','Plot','FontSize',fontSizeSmall,'Callback',{@plot_Callback2},'Value',0);
plotToggleTW2 = uicontrol('Parent',hPanel4,'Unit','Normalized','Position',[0.75 0.2 0.25 0.3],'Style','pushbutton','String','Clear','FontSize',fontSizeSmall,'Callback',{@cla_Callback2});

%%%%%%%%%%%%%%%%%%%%  Phase propagation plot panel %%%%%%%%%%%%%%%%%%%%%%%%
hPanel5 = uipanel('Title','Plot Phase Propagation','fontSize',fontSizeLarge,'Unit','Normalized','Position',[0.8 1-panelHeight 0.2 panelHeight]);

timeRangeprop0 = [-0.5 1];
uicontrol('Parent',hPanel5,'Unit','Normalized','Position',[0 0.5 0.5 0.5],'Style','text','String','Time Range (s)','FontSize',fontSizeMedium);
hTimeRangePropMin = uicontrol('Parent',hPanel5,'Unit','Normalized','BackgroundColor', backgroundColor,'Position',[0.5 0.5 0.25 0.5], ...
    'Style','edit','String',num2str(timeRangeprop0(1)),'FontSize',fontSizeSmall);
hTimeRangePropMax = uicontrol('Parent',hPanel5,'Unit','Normalized','BackgroundColor', backgroundColor,'Position',[0.75 0.5 0.25 0.5], ...
    'Style','edit','String',num2str(timeRangeprop0(2)),'FontSize',fontSizeSmall);

plotToggle = uicontrol('Parent',hPanel5,'Unit','Normalized','Position',[0 0 0.5 0.5],'Style','togglebutton','String','Plot/Pause','FontSize',fontSizeMedium,'Callback',{@plot_Callback2},'Value',0);
uicontrol('Parent',hPanel5,'Unit','Normalized','Position',[0.5 0 0.5 0.5],'Style','pushbutton','String','Clear','FontSize',fontSizeMedium,'Callback',{@cla_Callback2});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Handles for plot panels %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
numSelectedElectrodes = length(selectedElectrodes);
colorNamesElectrodes = jet(numSelectedElectrodes);

hGridPlots = getPlotHandles(1,2,[0.025 0.65 0.3 0.2],0.01,0.02,0); % for RF plots and electrode locations in MEA/ electrode locations and PSD's in EEG 
hTF = getPlotHandles(numSelectedElectrodes+1,2,[0.025 0.05 0.3 0.55],0.025,0.025); % for TF plots for all trials on the left and single trials on the right. Row1) All chan averaged. Row2-4) selected electrodes
hSignal = getPlotHandles(numSelectedElectrodes+1,numFrequencyRanges,[0.35 0.05 0.425 0.55],0.025,0.025);
hStats = getPlotHandles(1,numFrequencyRanges,[0.35 0.625 0.425 0.225],0.025,0.01);
hGridPlots2 = getPlotHandles(3,numFrequencyRanges,[0.8 0.05 0.19 0.8]);

%% start plotting and get additional required parameters
electrodeArray = showRFPositionsSelectedElectrodes(hGridPlots,goodElectrodes,selectedElectrodes,rfData,parameters,colorNamesElectrodes);

% plotting scripts 
%1. TF plots - done
% 2. PSD's -done
% 3. Bursts - done
% 4. PGD's - pending
% 5. Time series -done 
% 6. Polar plots - done
% 7. Gradient plots - pending
% 8. Phase propagation plots - done
% 9. RF positions - done
% 10. Topoplots - done
