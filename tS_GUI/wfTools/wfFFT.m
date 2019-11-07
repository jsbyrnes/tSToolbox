function [ F ] = wfFFT( T )
% this function takes an irisFetch trace structure (or array of structures)
% and returns a similar (array of)
% structure with added fields that contain the amplitude and phase spectra 
% of each trace as well as the frequencies corresponding to the points on 
% the spectra.
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

    FreqPerPoint = T(k).sampleRate/(T(k).sampleCount+1);
    
    pNyq=ceil(T(k).sampleCount/2);   % For even sampleCount this will correspond to Nyquist
                                           % For odd, it is 0.5 below Nyq (there is no Nyq)
    F(k).frequencies = (0:pNyq-1)*FreqPerPoint;    % freq values up and inc Nyq (in Hz)
    
    DFT = fft(detrend(T(k).data));
    
    ampSpec=abs(DFT);
    phSpec=angle(DFT);
    
    F(k).ampSpectrum=ampSpec(1:length(F(k).frequencies));
    F(k).phSpectrum=phSpec(1:length(F(k).frequencies));
    
end



