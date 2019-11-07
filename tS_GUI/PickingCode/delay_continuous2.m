function [ time_series ] = delay_continuous2( time_series, new_inc, delay)
%delay_continuous Applies a time shift of delay to time_series by up
%sampling, shifiting by a microsample, and downsampling again. 
    
    %upsample the time series
    
    %how much do you need to upsample by? Make it double to be safe
    upsample_rate = ceil(1/new_inc);%the factor of 2 is to ensure I am not aliasing. Not 100% sure it is needed.

    lag = round(delay*upsample_rate);
    
    time_series = resample(time_series, upsample_rate, 1);
    
    len = length(time_series);
    
    %apply the shift
    pad1=zeros(1,lag);
    time_series=[pad1 time_series ]; %makes the array at least as long as auxwf
    time_series=time_series(1:len); %shortens it if necessary
    
    %down sample to the original length
    
    time_series = resample(time_series, 1, upsample_rate);

end

