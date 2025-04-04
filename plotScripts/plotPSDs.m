function plotPSDs(data,timeVals,flag,taper,freqRangeList,selectedChannels,plotAxis)
% inputs 
% flag - 1 - get delta FFT
%        2 - uncorrected
% data - chan by trials by timevals
% taper - set number of tapers
% freqRangeList - cell of list of frequency ranges for visualization
% selectedChannels - list of required channels for plotting, set it as '0'
% to plot average across all electrodes
% plotAxis - plot handle

% Outputs 
% f - frequency values
% dFFT - power values
%% get fft
Fs = 1/(timeVals(2)-timeVals(1));
numChannels = size(data,1);
baseline = [-0.5 0];
stimPeriod = [0.25 0.75];
% fBL = 0:1/(diff(baseline)):Fs-1/(diff(baseline));
% fSt = 0:1/(diff(stimPeriod)):Fs-1/(diff(stimPeriod));
baseline = dsearchn(timeVals',baseline(1)):dsearchn(timeVals',baseline(2));
stimPeriod = dsearchn(timeVals',stimPeriod(1)):dsearchn(timeVals',stimPeriod(2));

params.tapers = [taper 1];
params.pad = -1;
params.Fs = Fs;
params.fpass = [0 200];
params.trialave = 1; % averaging across trials

for i = 1:numChannels
    % run fft 
    % logMeanFFTST = log10(mean(abs(fft(squeeze(data(i,:,baseline)),[],2))));
    % logMeanFFTBL = log10(mean(abs(fft(squeeze(data(i,:,stimPeriod)),[],2))));
    % logMeanFFTBL(1)=logMeanFFTBL(2);
    % dFFT(i,:) = logMeanFFTST - logMeanFFTBL;
    
    % run in Chronux
    [sStim,f] = mtspectrumc(squeeze(data(i,:,stimPeriod))',params);
    [sBase,~] = mtspectrumc(squeeze(data(i,:,baseline))',params);
    if flag==1
        dFFT(i,:) = log10(sStim)-log10(sBase);
    else
        dFFT(i,:) = log10(sStim);
    end
end

% plot the results per requirement
numFrequencyRanges = length(freqRangeList);
colorNames = gray(numFrequencyRanges);

if selectedChannels~=0
    colorVals = parula(numel(selectedChannels));
    for i = 1:numel(selectedChannels)
        plot(plotAxis,f,dFFT(selectedChannels(i)),'LineWidth',1.2,'Color',colorVals(i,:))
        hold(plotAxis,'on')
    end
else
    plot(plotAxis,f,mean(dFFT,1),'LineWidth',1.2,'Color','black')
    hold(plotAxis,'on')
end

for i = 1:numel(freqRangeList)
    xline(freqRangeList{i}(1),'Color',colorNames(i,:),'Parent',plotAxis)
    hold(plotAxis,'on')
    xline(freqRangeList{i}(2),'Color',colorNames(i,:),'Parent',plotAxis)
end


end

