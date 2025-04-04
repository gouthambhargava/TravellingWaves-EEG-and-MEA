% Modified from getBurstLengthHilbert in https://github.com/supratimray/GammaLengthProjectCodes.

function [burstTS] = getHilbertBurstER(bpfSignal,timeVals,thresholdFactor,stimulusPeriodS,baselinePeriodS,analysisPeriodS)

if ~exist('stimulusPeriodS','var');     stimulusPeriodS=[0.25 0.75];    end
if ~exist('baselinePeriodS','var');     baselinePeriodS=[-1 0];         end
if ~exist('analysisPeriodS','var');     analysisPeriodS=stimulusPeriodS;end

hilbertPower=abs(hilbert(bpfSignal'))';

mBL = getMedianBaseline(hilbertPower,timeVals,baselinePeriodS);
% if isempty(thresholdFactor)
%     mST = getMedianBaseline(hilbertPower,timeVals,stimulusPeriodS);
%     thresholdFactor = mean(mST./mBL);
%     disp(['Using threshold factor of: ' num2str(thresholdFactor)]);
% end

threshold=thresholdFactor*mBL;
hilbertPower = reshape(hilbertPower,[1,size(hilbertPower,1),size(hilbertPower,2)]);
numTrials = size(hilbertPower,2);
for i=1:numTrials
    burstTS = getBurstLengthHilbertSingleTrial(hilbertPower(:,i,:),timeVals,threshold,analysisPeriodS);
end
end

function [burstTS] = getBurstLengthHilbertSingleTrial(bandPassPowerSingleTrial,timeVals,threshold,analysisPeriodS)

stPos = intersect(find(timeVals>=analysisPeriodS(1)),find(timeVals<analysisPeriodS(2)));
st=timeVals(1,stPos);
stBandPassPowerSingleTrial=bandPassPowerSingleTrial(1,stPos);

seedPosList=find((stBandPassPowerSingleTrial/threshold)>1); % Seed point is any time point which exceeds the threshold
burstStartPosList = [];
burstEndPosList = [];
timePos=1;
timeLength=length(st);
while (timePos < timeLength) % As long as the time position is less than the search interval
    tPos = seedPosList(find(seedPosList>timePos,1)); % Find the position of the next seed
    minPos = find((stBandPassPowerSingleTrial(1:tPos)/(threshold/2))<1,1,'last');
    if ~isempty(tPos)
        maxPos = find((stBandPassPowerSingleTrial(tPos+1:timeLength)/(threshold/2))<1,1,'first');

        if isempty(minPos)
            tMinPos=1;
        else
            tMinPos=minPos;
        end
        if isempty(maxPos)
            tMaxPos=timeLength;
        else
            tMaxPos=tPos+maxPos;
        end
        timePos=tMaxPos;
        burstStartPosList=cat(1,burstStartPosList,stPos(tMinPos));
        burstEndPosList=cat(1,burstEndPosList,stPos(tMaxPos));
    else
        timePos=timeLength;
    end
end
% 
burstTS = nan(1,length(bandPassPowerSingleTrial));
for ind = 1:length(burstStartPosList)
    burstTS(burstStartPosList(ind):burstEndPosList(ind)) = 1;
end
end

function mBL = getMedianBaseline(smoothedPowerBpfSignal,timeVals,baselinePeriodS)
blPos = intersect(find(timeVals>=baselinePeriodS(1)),find(timeVals<baselinePeriodS(2)));
blData = smoothedPowerBpfSignal(:,blPos);
mBL = median(blData(:));
end
