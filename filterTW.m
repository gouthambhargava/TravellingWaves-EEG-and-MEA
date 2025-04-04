function [FS,filteredSignal,bandPhase, bandPower] = filterTW(data,srate,freqBands,analysisWindow,freqReq,timeVals)

% Instantaneous frequency estimates obtained using "Frequency Sliding" methods 
% http://mikexcohen.com/data/Cohen2014_freqslide.pdf

%Inputs:
% data - signal to analyze.  chan x time.
% srate - sampling rate of the signal in Hz
% freqbands - the required band for the analysis. Ex. [8 12];
% freqReq - if 1 - complete the frequency sliding script - uses equiripple
%                  filter (one way)
%              2 - Do hilbert burst detection - uses butterworth filter
%              (filtfilt)
%              3 - Return FS as a singleton matrix where each time point is considered as burst point      

% Outputs 
% frequency sliding- instantaneous frequency of signal in each band.
% filteredSignal - filtered signal in the specified band. 
% hilbertPhase - instantaneous phase of signal in each channel
% hilbertPower - instantaneous power of signal in each channel


%Key Steps.  Implementations of these steps are commented throughout code.

%    Calculates "frequency sliding" (MX Cohen 2014; JNeuroscience) in each
%    band.  Note that these frequency estimates are continuous values within a
%    band and are NOT rounded to wavefreqs.

%    Removes frequency estimates outside of detected band occuring due to
%    phase slips, see Cohen 2014 JNeurosci. Figure 1B and caption.

%% set up some parameters
% demean data
data= data-mean(data,2);
numElectrodes = size(data,1);

%% set filter parameters
transWidth    = .15; %width of transition zone
idealresponse  = [ 0 0 1 1 0 0 ];
filtfreqbounds = [ 0 (1-transWidth)*freqBands(1) freqBands(1) freqBands(2) freqBands(2)*(1+transWidth) srate/2 ]/(srate/2);
filtOrder     = round(2*(srate/freqBands(1))); 
filterweights  = firls(filtOrder,filtfreqbounds,idealresponse);

% initialize outputs
filteredSignal = nan(size(data,1),size(data,2)); 
bandPhase = nan(size(data,1),size(data,2)); %Initialize phases
bandPower = nan(size(data,1),size(data,2)); %Initialize phases

%% on to the main loop
for ind = 1:size(data,1) % loop across all channels
    signal = data(ind,:);
    % Frequency sliding code from MX Cohen.
    % filter data    
    filteredSignal(ind,:) = filtfilt(filterweights,1,signal);
    
    %hilbert the filtered signal
    hilbertPhase = angle(hilbert(filteredSignal(ind,:)));
    bandPhase(ind,:) = hilbertPhase;
    bandPower(ind,:) = abs(hilbert(filteredSignal(ind,:))).^2;
end

if freqReq ==1   % get frequency sliding estimates
        FS = nan(size(data,1),size(data,2)); %either frequency sliding estimates or filtered signal
        for ind = 1:size(filteredSignal,1) % loop across all channels
            %code from MX Cohen
            signal = filteredSignal(ind,:);
            instFreq = nan(1,length(signal));
            instFreq(1:end-1) = srate*diff(unwrap(bandPhase(ind,:)))/(2*pi); %code from fs paper cohen 2014
            time_wins = [.05 .2 .4]; %time windows in fractions of a second from MX Cohen
            orders = time_wins*srate;

            %window signal into 10 epochs to make it more tractable.
            numchunks = 10;
            chunks = floor(linspace(1,length(instFreq),numchunks)); %make epochs
    
            meds = zeros(length(orders),length(instFreq));
            for iWin = 1:length(orders)%median filter using different window sizes. 
                for iChunk = 2:numchunks
                    chunkidx = chunks(iChunk-1):chunks(iChunk)-1; %don't double count edges, last sample will be excluded.
                    meds(iWin,chunkidx) = medfilt1(instFreq(chunkidx),round(orders(iWin)));
                end
            end
    
            %take the median value across different medians
            median_of_meds = median(meds);
    
            % Key Step #4. NaN out frequency estimates outside of the filter band
            % clear below* above* outside*
            below_idx = (median_of_meds<freqBands(1));
            above_idx = (median_of_meds>freqBands(2));
            outside_idx = find(below_idx+above_idx==1);
            median_of_meds(outside_idx)=NaN;
            FS(ind,:) = median_of_meds; %all frequency sliding estimates within band.
        end
elseif  freqReq ==2 %get hilbert burst estiamtes
        % repeat the data matrix to create pseudo trials 
        FS = nan(size(data,1),size(data,2)); %either frequency sliding estimates or filtered signal

        if freqBands(1)<20
            error('Hilbert burst detection cannot be used for bands slower than gamma')
        end     
        thresholdFactor = 4;
        baselinePeriodS = [-0.5 0];
        stimulusPeriodS = [0.25 0.75];
        analysisPeriodS = [-0.5 1];
        filterOrder = 4;
        % for iElec=1:numElectrodes
        %     [~,~,~,FS(iElec,:,:)] = getHilbertBurst(squeeze(data(:,iElec,:)),timeVals,thresholdFactor,0,stimulusPeriodS,baselinePeriodS,freqBands,filterOrder,1,analysisPeriodS);
        % end
        [~,~,~,FS] = getHilbertBurst(data,timeVals,thresholdFactor,0,stimulusPeriodS,baselinePeriodS,freqBands,filterOrder,1,analysisPeriodS);

        % [FS] = getHilbertBurstER(data,timeVals,thresholdFactor,stimulusPeriodS,baselinePeriodS,analysisPeriodS);
else 
        FS = ones(size(data)); 
end
FS(~isnan(FS)) = 1;
analysisWindowIndices = dsearchn(timeVals',analysisWindow(1)):dsearchn(timeVals',analysisWindow(2));
FS(:,setdiff(1:length(FS),analysisWindowIndices)) = nan;
end
