function [ outTr ] = wfButterworth( inTr, freqs )
% wfButterworth applies a cosine-tapered filter to a waveform (IRISfetch trace) 
% structure and returns a similar structure with the data field now 
% containing the filtered trace
% USAGE:
% [ outTr ] = wfButterworth( inTr, freqs  )
% where freqs is the 1x2 vector [L1 H1]
% L1 is low-frequency side of low-frequency taper
% H1 is low-frequency side of  hi-frequency taper
% all of the above should be specified in Hz

%the following is just a fancy way of getting the elements of freqs into 4
%variables.
fCell=num2cell(freqs);
[L1, H1] = fCell{:};

outTr=inTr;

%d = designfilt('bandpassiir', 'FilterOrder', 12, 'HalfPowerFrequency1', ...
%     L1, 'HalfPowerFrequency2', H1, 'SampleRate', ...
%     inTr(1).sampleRate, 'DesignMethod', 'butter');

d = designfilt('bandpassiir', ...
    'FilterOrder', 20, 'PassbandFrequency1', L1, 'PassbandFrequency2', H1,...
    'StopbandAttenuation1', 30, 'PassbandRipple', 0.1, 'StopbandAttenuation2', 30, 'SampleRate', inTr(1).sampleRate);

for i=1:length(inTr)
    
    ts = inTr(i).data;
    n = length(ts);
    
    tap = tukeywin(n, 1/(n*L1/inTr(i).sampleRate));
    
    %minimum phase filtering    
    %outTr(i).data = filter(d, ts.*tap);
    
    %zero phase filtering    
    outTr(i).data = filtfilt(d, ts.*tap);
    
end

    %or zero phase filtering, simple
    %[z,p,k] = butter(3,[L1 H1].*(2/inTr(i).sampleRate));
    %[z,p,k] = ellip(3, 0.001, 100, [low high].*(2*dt));
    %[z,p,k] = cheby2(3, 50, [low high].*(2*dt));
    %[SOS, G] = zp2sos(z,p,k);
            
    %Apply the filter to the time series and we're done
    %outTr(i).data = filtfilt(SOS, G, ts.*tap);

