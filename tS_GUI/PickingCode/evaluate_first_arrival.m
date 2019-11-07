function [ ind, s2n, s2n_secondary, trigger, secondary_trigger] = evaluate_first_arrival( trace, sta_window, lta_window, trigger_level, search_ind, secondary_trace, secondary_trigger_level)
%evaluate_arrival Grid searchs over the trace to find where the trace is
%best split into high/low amplitude

    if nargin == 4
        
        search_ind = lta_window + 1:(length(trace) - sta_window - 1);
               
    end

    trigger = 0;
    secondary_trigger = 0;
    s2n = 0;
    
    if nargin ~= 7
    
        secondary_trigger_level = 0;
        
    end
    
    %do a lta/sta picker    
    for i = search_ind%don't check the extreme end points
       
        lta = trace(i - lta_window:i);
        
        if (i+sta_window) >= length(trace)
            
            break
            
        end
        
        sta = trace(i:i+sta_window);
        
        lta = rms(lta);
        sta = rms(sta);
               
        trigger(i) = sta/lta;
        
        if nargin == 7
            
            lta = secondary_trace(i - lta_window:i);
            sta = secondary_trace(i:i+sta_window);
            lta = rms(lta);
            sta = rms(sta);
               
            secondary_trigger(i) = sta/lta;
            
        end
        
        if trigger(i) >= max(trigger)
            
            ind = i;
            s2n = trigger(i);
            
        end
                    
    end
    
    %smooth the trigger over the STA
    trigger = smooth(trigger, sta_window);
    secondary_trigger = smooth(secondary_trigger, sta_window);
    
    if max(secondary_trigger) < secondary_trigger_level || max(trigger) < trigger_level
    
        ind = [];
        s2n = [];
        s2n_secondary = [];
        
    else
        
        s2n_secondary = max(secondary_trigger);
        
    end
    
end

