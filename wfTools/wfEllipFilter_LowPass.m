function [ outTr, d ] = wfEllipFilter_LowPass( inTr, freqs, phase, d)
% wfButterworth applies a cosine-tapered filter to a waveform (IRISfetch trace) 
% structure and returns a similar structure with the data field now 
% containing the filtered trace
% USAGE:
% [ outTr ] = wfButterworth( inTr, freqs  )
% where freqs is the 1x2 vector [L1 H1]
% L1 is low-frequency side of low-frequency taper
% H1 is low-frequency side of  hi-frequency taper
% all of the above should be specified in Hz

outTr=inTr;

if nargin ~= 4

d = designfilt('lowpassiir', 'PassbandFrequency', freqs-0.05, 'StopbandFrequency', freqs+0.05, ...
    'PassbandRipple', 1, 'StopbandAttenuation', 60, 'SampleRate', inTr(1).sampleRate, 'DesignMethod', 'ellip');

end

for i=1:length(inTr)
    
    ts = inTr(i).data;
    n = length(ts);
    
    tap = tukeywin(n, 1/(n*freqs(1)/inTr(i).sampleRate));
    
    if nargin < 3
    
        %zero phase filtering    
        outTr(i).data = filtfilt(d, ts.*tap);
    
    elseif nargin >= 3
    
        if strcmp(phase, 'zero')
        
            %zero phase filtering    
            outTr(i).data = filtfilt(d, ts.*tap);
            
        elseif strcmp(phase, 'minimum')
           
            %minimum phase filtering    
            outTr(i).data = filter(d, ts.*tap);
            
        else
            
            error('Invalid phase flag');
            
        end
    
    end
    
end

    %or zero phase filtering, simple
    %[z,p,k] = butter(3,[L1 H1].*(2/inTr(i).sampleRate));
    %[z,p,k] = ellip(3, 0.001, 100, [low high].*(2*dt));
    %[z,p,k] = cheby2(3, 50, [low high].*(2*dt));
    %[SOS, G] = zp2sos(z,p,k);
            
    %Apply the filter to the time series and we're done
    %outTr(i).data = filtfilt(SOS, G, ts.*tap);

