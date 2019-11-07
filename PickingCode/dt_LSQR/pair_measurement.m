function [ dt, flip ] = pair_measurement( To, Ts )
%gives the pairwise t* for two traces

    %try measuring this way
    flip = 0;
    
    To.data = To.data/max(To.data);
    Ts.data = Ts.data/max(Ts.data);
       
    [cc, lags]   = xcorr(To.data, Ts.data);
    [~, tmpind]  = max(cc);
    dt           = lags(tmpind)/To.sampleRate;
            
end

