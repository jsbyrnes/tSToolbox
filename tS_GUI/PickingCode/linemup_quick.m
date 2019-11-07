function [ swf ] = linemup_quick( swf, owf )
% swf=synthetic waveform, owf=observed waveform

    delay = delayest_3point(owf,swf);
        
    [ swf ] = delay_continuous( swf', 1, delay )';

end