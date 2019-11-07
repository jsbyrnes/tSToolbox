function [ outTr ] = wfCosFilter( inTr, freqs  )
% wfLoPass applies a cosine-tapered filter to a waveform (IRISfetch trace) 
% structure and returns a similar structure with the data field now 
% containing the filtered trace
% USAGE:
% [ outTr ] = wfCosFilter( inTr, freqs  )
% where freqs is the 1x4 vector [L1 L2 H1 H2]
% L1 is low-frequency side of low-frequency taper
% L2 is  hi-frequency side of low-frequency taper
% H1 is low-frequency side of  hi-frequency taper
% H2 is  hi-frequency side of  hi-frequency taper
% all of the above should be specified in Hz

%the following is just a fancy way of getting the elements of freqs into 4
%variables.
fCell=num2cell(freqs);
[L1, L2, H1, H2] = fCell{:};

outTr=inTr;

for k=1:length(inTr)
    
    Nfft = length(inTr(k).data);         % number  of points in fft = n of samples
    dF = inTr(k).sampleRate/Nfft;       % frequency interval
    
    Nyq=(inTr(k).sampleRate/2);        % Nyquist frequency
                                           
    freqVals = (0:(Nfft-1))*dF;         % this gives us frequencies with the correct spacing but going all the way to the sampling frequency
    
    %this makes the negative frequencies such that they are where they 
    %should be with respect to the complex amplitudes in F.
    freqVals(freqVals>Nyq) = freqVals(freqVals>Nyq)-(Nyq*2);        

    %to make specifying the pass band easier, I make all frequencies
    %postive
    freqVals=abs(freqVals);


    %A is the amplitude spectrum of the filter, we start by setting it all
    %to 1, only values between L2 and H1 will be 1 at the end
    A = ones(length(inTr(k).data),1);
    %outside range
    A(freqVals<L1 | freqVals>H2) = 0;
    %up taper
    pnt = find(freqVals>=L1 & freqVals<=L2);
    A(pnt) = 0.5*(1-cos(pi*(freqVals(pnt)-L1)/(L2-L1)));
    %down taper
    pnt = find(freqVals>=H1 & freqVals<=H2);
    A(pnt) = 0.5*(1+cos(pi*(freqVals(pnt)-H1)/(H2-H1)));

    
    %make sure the trace starts at about zero
    d=inTr(k).data;
    inTr(k).data=inTr(k).data-mean(d(1:round(inTr(k).sampleRate)));% take out the mean of the first second
    
    %Apply the filter to the time series and we're done
    outTr(k).data=real(ifft(fft(inTr(k).data).*A));
    

end

