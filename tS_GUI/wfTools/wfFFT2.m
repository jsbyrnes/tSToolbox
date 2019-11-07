function [ F ] = wfFFT2( T )
% this function takes an irisFetch trace structure (or array of structures)
% and returns a similar (array of)
% structure with added fields that contain the amplitude and phase spectra 
% of each trace as well as the frequencies corresponding to the points on 
% the spectra. This function returns the two-sided (positive and neg freqs)
% spectra.
%
%
% USAGE 
% [ F ] = wfFFT( T )
% where T is a trace structure.
% the new fileds in F will be ampSpectrum, phSpectrum, frequencies.
% to plot the amplitude spectrum do:
% plot(F(k).frequencies,F(k).ampSpectrum)
% where k is an integer greater than zero and smalle than the length of the
% 

F=T;

for k=1:length(T)

    
    if T(k).sampleCount ~= length(T(k).data)
        %disp('Warning: sampleCount did not match length of data, was reset')
        T(k).sampleCount=length(T(k).data);
    end
    
    dFreq = T(k).sampleRate/T(k).sampleCount; %frequency spacing
    
    fNyq=T(k).sampleRate/2;   %Nyquist frequency
    
    %the next two steps build the frequency vector (freqs at which spectrum
    %was calculated) postivie and negative
    fVector=(0:T(k).sampleCount-1)*dFreq;
    fVector(fVector>fNyq)=fVector(fVector>fNyq)-fNyq*2;
    
    F(k).frequencies = fVector; 
    
    DFT = fft(detrend(T(k).data));
    
    ampSpec=abs(DFT);
    ampSpec=ampSpec/T(k).sampleCount;
    phSpec=angle(DFT);
    
    F(k).ampSpectrum=ampSpec(1:length(F(k).frequencies));
    F(k).phSpectrum=phSpec(1:length(F(k).frequencies));
    
end



