function [ al_tr ] = linemup_fast( swf, owf )
% swf=synthetic waveform, owf=observed waveform
    
    madeSwitch=false;
    %make rows if colummns
    if size(swf,1) == length(swf)
        madeSwitch=true;
        swf=swf';
        owf=owf';
    end

    owf = [ zeros(size(owf)) owf ];%shift forward
    swf = [ swf zeros(size(swf)) ];
    
    [al_tr, ~] = alignsignals(swf, owf, length(owf), 'truncate');
        
    al_tr = al_tr(round(length(al_tr)/2) + 1:end);
    
    if madeSwitch;
       al_tr=al_tr'; 
    end

end
