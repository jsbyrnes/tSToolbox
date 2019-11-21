function [Aratio, SNR] = calcScatterIndex(traces,filterFreqs, arrival, sec_before, sec_after, before_start, after_start)
%the idea here is to take the envelope of the trace (tapered) and a
%high-pass filtered trace (tapered) and take the ratio of the area under
%each.
% USAGE: 
% [Aratio] = calcScatterIndex(traces,filterFreqs)
% the output is a vector, suitable for adding to colorby
% traces input should be a traces structure
% filterFreqs has 2 corner frequencies for a butterworth bandpass filter

for k = 1:length(traces)

    %trace and filtered trace
    tr          = traces(k);
    arrival_ind = arrival*tr.sampleRate;
    
    tr.data = detrend(tr.data); %remove mean of first second
    tr.data = tr.data - mean(tr.data(1:tr.sampleRate));
    %trF=wfCosFilter(tr,filterFreqs);
    trF     = wfButterworth(tr,filterFreqs);

    %now build a taper and apply it, 1 second on each side, linear rolloff for
    %now.
    rollOff                    = tr.sampleRate; %how many samples in 1 second
    taper                      = ones(size(tr.data));
    taper(1:rollOff)           = linspace(0,1,rollOff);
    taper(end-(rollOff-1):end) = linspace(1,0,rollOff);

    tr.data  = tr.data.*taper;
    trF.data = trF.data.*taper;

    %calculate envelopes
    trEnv  = abs(hilbert(tr.data));
    trFEnv = abs(hilbert(trF.data));

    before = (arrival_ind - (before_start + sec_before)*tr.sampleRate:arrival_ind - before_start*tr.sampleRate); % 5 before the arrival, take a second out
    after  = (arrival_ind + tr.sampleRate*after_start:arrival_ind + tr.sampleRate*(after_start + sec_after)); % ten seconds after the arrival

    A      = trapz(trEnv(after)); % Area under the envelope using the trapezoid method
    AF     = trapz(trFEnv(after)); % Area under the envelope using the trapezoid method
    Anoise = trapz(trFEnv(before)); % Area under the envelope of the noise
    AF     = AF-(sec_after/sec_before)*Anoise; % the 2 is beacause I'm comparing 5 s to 10 s.
    AF     = max(AF,0);

    Aratio(k) = AF/A;
    SNR(k)    = AF/Anoise;
        
end

