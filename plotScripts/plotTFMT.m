function plotTFMT(axis,data,timeVals,axisRanges,freqRangeHz,colorNames,type,trialNo,stimulusPeriodS)
% supply data in the form channels x trials x times
% set trialNo to [] if all channels are needed

Fs = round(1/(timeVals(2)-timeVals(1)));

if ~isempty(trialNo)
    data = data(:,trialNo,:);
end

% TF options
movingwin = [0.25 0.025];
params.tapers = [1 1];
params.pad = -1;
params.Fs = Fs;
params.fpass = [0 200];
params.trialave = 1; % averaging across trials
blRange = [-0.5 0];

% Time frequency analysis with multitapers
for j = 1:size(data,1)
    [S,timeTF,freqVals] = mtspecgramc(squeeze(data(j,:,:))',movingwin,params);
    logPower(:,:,j) = log10(S);
end
xValToPlot = timeTF+timeVals(1)-1/Fs;

logPower = mean(logPower,3);

if strcmp(type,'delta')
    blPos = intersect(find(xValToPlot>=blRange(1)),find(xValToPlot<blRange(2)));
    logBLPower = repmat(mean(logPower(blPos,:,:),1),length(xValToPlot),1);
    deltaPower = 10*(logPower - logBLPower);
    pcolor(hTF,xValToPlot,freqVals,deltaPower');
else
    pcolor(hTF,xValToPlot,freqVals,logPower');
end

shading(hTF,'interp');
colormap(hTF,'jet');
axis(hTF,[axisRanges{2} axisRanges{1}]);
clim(hTF,axisRanges{3});

% Indicate frequency ranges of interest
numRanges = length(freqRangeHz);
for j=1:numRanges
    line([axisRanges{2}],[freqRangeHz{j}(1) freqRangeHz{j}(1)],'color',colorNames(j,:),'parent',hTF);
    line([axisRanges{2}],[freqRangeHz{j}(2) freqRangeHz{j}(2)],'color',colorNames(j,:),'parent',hTF);
end

if ~isempty(stimulusPeriodS)
    line([stimulusPeriodS(1) stimulusPeriodS(1)],[axisRanges{1}],'color','k','parent',hTF);
    line([stimulusPeriodS(2) stimulusPeriodS(2)],[axisRanges{1}],'color','k','parent',hTF);
end
end