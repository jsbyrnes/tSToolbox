function [ tS, f, R, s2n1, s2n2] = spectral_slopes_tS( fs, f_range, t1, t2, flag, Ed )
%spectral_slopes_tS t is the time axis, t1 and t2 are two time series, 
%flag indicates how the smoothing is handled, and 1 variable argument can
%be used for certain types of smoothing. tS out is the attenuation of t1
%relative to t2
    
    f = (fs*(0:(length(t1)/2))/length(t1))';

    if strcmp(flag, 'MT')

        A1 = zeros(size(t1));
        A2 = zeros(size(t2));

        [~,K] = size(Ed);
        
        for i = 1:K

            A1 = A1 + fft(t1.*Ed(:,i));
            A2 = A2 + fft(t2.*Ed(:,i));

        end

        A1 = abs(A1/K);
        A2 = abs(A2/K);
        
    else

        A1 = abs(fft(t1));
        A2 = abs(fft(t2));

    end
        
    A1 = A1(1:floor(length(t1)/2) + 1);
    A2 = A2(1:floor(length(t2)/2) + 1);
            
    f_ind = (f >= f_range(1)) & (f <= f_range(2));
    f = f(f_ind);
    
    %get the ln ratio of the amplitudes
    R = log(A1(f_ind)./A2(f_ind));
    %get the slope and solve for tS (t*)
    p = polyfit(f, R, 1);
    tS = -p(1)/pi;
    
end

