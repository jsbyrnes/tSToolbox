function [ al_tr, lag_index, minor_lag_index ] = linemup_jsb( swf, owf, fitting_window )
% swf=synthetic waveform, owf=observed waveform

    fine_search_inc   = 0.0001;%for the deriviatives in the finer search gridsearch, and the upsamplin rate for delaycontinuous2
    fine_search_start = 0.5;%half width of the finer search
    terminate_at      = 0.0000001;%level to terminate at
    fine_search_inter = 10000;%how many euler's method steps to take
    
    madeSwitch=false;
    %make rows if colummns
    if size(swf,1) == length(swf)
        madeSwitch=true;
        swf=swf';
        owf=owf';
    end

    swf(1:100)=[]; % this gets rid of the first 100 samples. You are only
    %testing 'forward' lags

    %jsb commented out 
%     maxlag=200;
%     linerr=nan(1,maxlag);
%     atrM=nan(maxlag,length(owf));

    lagvec = 0:400;

    for lag = 1:length(lagvec)
        
        pad1=zeros(1,lag);
        pad2=zeros(1,length(owf)-(length(swf)+lag));
        a1=[pad1 swf pad2]; %makes the array at least as long as auxwf
        a1=a1(1:length(owf)); %shortens it if necessary

        %JSB edit
        %use delayseq instead of shifting by samples to align the time series
        %account for time shifts of less than a sample
        %The input variable is a1 if you pad and swf if not
        %a1 = delay_continuous(a1, 1, lagvec(lag));%just do it in samples, so set fs=1

        atrM(lag,:)=a1;

        try
        
            %linerr(lag)=norm(a1-owf);
            linerr(lag)=norm(a1(fitting_window(1):fitting_window(2)) - owf(fitting_window(1):fitting_window(2)));
        
        catch
            
            keyboard
        
        end
        
    end

    [~, shft1]=min(linerr);
    al_tr=atrM(shft1,:);

    lag_index = lagvec(shft1);
    minor_lag_index = 0;
%     xc = xcorr(owf, swf);
%     [~, ind] = max(xc);
%     lag_index = ind - length(owf);
% 
%     pad1=zeros(1,lag_index);
%     pad2=zeros(1,length(owf)-(length(swf)+lag_index));
%     a1=[pad1 swf pad2]; %makes the array at least as long as auxwf
%     al_tr=a1(1:length(owf)); %shortens it if necessary
    
    fine_shift = fine_search_start;%staring value
    
    %al_tr = [ al_tr(fine_search_width + 1: end) zeros(1, fine_search_width) ];
    
    fitting_window = fitting_window + lag_index;%move it over so it still works
    
    for i = 1:fine_search_inter
       
        %get the deritvate
                
        a1 = delay_continuous(al_tr, 1, fine_shift);%just do it in samples, so set fs=1
        a2 = delay_continuous(al_tr, 1, fine_shift + fine_search_inc);%just do it in samples, so set fs=1
        a3 = delay_continuous(al_tr, 1, fine_shift - fine_search_inc);%just do it in samples, so set fs=1
%         a1 = delay_continuous2(al_tr, 20, fine_shift);%just do it in samples, so set fs=1
%         a2 = delay_continuous2(al_tr, 20, fine_shift + fine_search_inc);%just do it in samples, so set fs=1
%         a3 = delay_continuous2(al_tr, 20, fine_shift - fine_search_inc);%just do it in samples, so set fs=1

        %get the central difference
        linerr1=norm(a1 - owf);
        linerr2=norm(a2 - owf);
        linerr3=norm(a3 - owf);
        
        first_deriv = (linerr2 - linerr3)/(2*fine_search_inc);
        
        %for Euler's method
        %fine_shift = fine_shift - fine_search_step*first_deriv;
        
        %for newton's method
        secderiv    = (linerr2 - 2*linerr1 + linerr3)/fine_search_inc^2;
        
        if secderiv ~= 0
        
            step_size = -1*first_deriv/secderiv;
            
        else
            
            step_size = 0; %too close, code will give nans and spin
            
        end
        
        if abs(first_deriv/secderiv) > 0.1
            
            %disp('Warning - step size in refinement too large, probably linear misfit');
            %will do Euler's method instead
            
            step_size = -(first_deriv)*0.1;
            
        end
                
        fine_shift = fine_shift + step_size;
                
        %you can see how it changed over time
        shift_history(i) = fine_shift;
        
        if i > 1 && abs(shift_history(i - 1) - shift_history(i)) < terminate_at
            
            break%close enough
            
        end
                
    end

%     if i == fine_search_inter
%         
%         disp('Warning - lineemup not converging');
%         
%     end
    
%     al_tr = delay_continuous2(al_tr, 20, fine_shift);
    al_tr = delay_continuous(al_tr, 1, fine_shift);
    
    minor_lag_index = fine_shift - fine_search_start/2;
    
    %you can probably do this without a for loop 

    if madeSwitch;
       al_tr=al_tr'; 
    end

end

%%%%%
%old code grid searches over 
    %exact same code, but redo within 1 sample to get a better measurement
    
    %need to shift back by one sample to use these lags
%     lagvec = 0:fine_search_inc:2*fine_search_width;%small, subsampling frequency lags
%     
%     linerr = [];
%     atrM   = [];
%     
%     al_trtmp = [ al_tr(fine_search_width + 1: end) zeros(1, fine_search_width) ];
%     
%     for lag2 = 1:length(lagvec)
% 
%         a1 = delay_continuous2(al_trtmp, fine_search_inc, lagvec(lag2));%just do it in samples, so set fs=1
% 
%         atrM(lag2,:)=a1;
% 
%         try
%             
%             %see comments from above
%             %linerr(lag2)=norm(a1-owf);
%             linerr(lag2)=norm(a1(fitting_window(1):fitting_window(2)) ...
%                 - owf(fitting_window(1):fitting_window(2)));
%         
%         catch
%             
%             keyboard
%         
%         end
%         
%     end
% 
%     [m, shft2]=min(linerr);
%     al_tr=atrM(shft2,:);


%%%%
%attempt to find the minimum with Newton's method - but the second
%derivative is too small until you get kind close
%     fine_shift = fine_search_width;%staring value
%     
%     for i = 1:fine_search_inter
%        
%         %get the deritvate
%         a1 = delay_continuous2(al_trtmp, fine_search_inc, fine_shift);%just do it in samples, so set fs=1
%         a2 = delay_continuous2(al_trtmp, fine_search_inc, fine_shift + fine_search_inc);%just do it in samples, so set fs=1
%         a3 = delay_continuous2(al_trtmp, fine_search_inc, fine_shift - fine_search_inc);%just do it in samples, so set fs=1
% 
%         %get the central difference
%         linerr1=norm(a1(fitting_window(1):fitting_window(2)) ...
%             - owf(fitting_window(1):fitting_window(2)));
%         linerr2=norm(a2(fitting_window(1):fitting_window(2)) ...
%             - owf(fitting_window(1):fitting_window(2)));
%         linerr3=norm(a3(fitting_window(1):fitting_window(2)) ...
%             - owf(fitting_window(1):fitting_window(2)));
%         
%         first_deriv = (linerr2 - linerr3)/(2*fine_search_inc);
%         secderiv    = (linerr2 - 2*linerr1 + linerr3)/fine_search_inc^2;
%         fine_shift = fine_shift - first_deriv/secderiv;
%         
%         shift_history(i) = fine_shift;
%         
%         %al_tr = delay_continuous2(al_tr, fine_search_inc, fine_shift);
%         
%     end


