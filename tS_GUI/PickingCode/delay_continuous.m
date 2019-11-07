function [ time_series ] = delay_continuous( time_series, fs, delay )
%delay_continuous Applies a time shift of delay to time_series in the
%frequency domain. Almost works perfectly but there is a small imaginary
%component coming out so I am doing something wrong, most likely related to
%how fftshift is being used. 
    
    undo = 0;

    %you need time_series to be odd
%     if mod(length(time_series), 2) == 0
%     
%         time_series(end + 1) = 0;
%         undo                 = 1;
%         
%     end
        
    dFreq = fs/length(time_series); %frequency spacing
    
    fNyq=fs/2;   %Nyquist frequency
    
    %the next two steps build the frequency vector (freqs at which spectrum
    %was calculated) postivie and negative
    f=(0:length(time_series)-1)*dFreq;
    f(f>fNyq)=f(f>fNyq)-fNyq*2;
    
    %f = linspace(-1*fs/2, fs/2, length(time_series));
    
    Ftime_series = fft(time_series);

    Ftime_series = Ftime_series .* exp (-1 * 1i * 2 * pi * f * delay);
    
    time_series  = ifft((Ftime_series), 'symmetric');
    
    if undo
        
        time_series = time_series(1:end-1);
        
    end
    
end

