function [outTraces] = wfResample(inTraces,newSR)
%USAGE: [outTraces] = wfResample(inTraces,newSR);
% inTraces, outTraces are irisFetch type structures of arrays
% newSR is the desired sampling rate in sps.
% in outTraces all the .data will have the same number of samples
% the fields sampleCount and sampleRate will be updated accordingly.
% I should add a comment saying the traces were resampled.
% any traces that are shorter than what they are supposed to be, will be
% axed.

%get rid of 'too short' traces
ttime=[inTraces.endTime]-[inTraces.startTime]; %length of each trace
ttime=ttime*24*60*60; %in seconds now
tDiff=max(ttime)-ttime; %max ttime is the nominal time, there's no reason to have any that are too long

oneSample=1./[inTraces.sampleRate];
ix=tDiff > oneSample*1.5; %find the ones that are off by more than one and a half samples

%get rid of them
inTraces=inTraces(~ix);


%make sure newSR is an integer
if ~(newSR==round(newSR))
    error('newSR needs to be an integer')
end

%loop through traces and resample as needed


outTraces=inTraces;
for k=1:length(inTraces)
    oldSR=inTraces(k).sampleRate;
    
    %check if resampling is needed
    if ~(oldSR==newSR)
        
        %check if oldSR is an integer, otherwise make it one and adjust
        %newSR accordingly. Resample
        if (oldSR==round(oldSR))
            %it's an integer
            outTraces(k).data=resample(inTraces(k).data,newSR,oldSR);

        else
            %it's not an integer
            oldSR1=round(oldSR*1000); %keep 3 decimals
            newSR1=newSR*1000;
            outTraces(k).data=resample(inTraces(k).data,newSR1,oldSR1);
            
        end
        
        outTraces(k).sampleRate=newSR;
        outTraces(k).sampleCount=length(outTraces(k).data);
                
    end
    
    outTraces(k).data = outTraces(k).data - mean(outTraces(k).data);
    
end

%OK, now go back and make sure that the sampleCounts are all the same.
usc=unique([outTraces.sampleCount]);

if length(usc)==1
    %if they are already uniform
    return
elseif length(usc)>2 || (max(usc)-min(usc)) > 1
    %the expected scenario is that there will be no more than two sample
    %counts and they will be X and X+1, if that's not the case, display
    %this message
    warning('something is not as expected, check things carefully')
    keyboard
end

%if you've made it this far then there's more than one sampleCount
%go through and fix that
minSc=min(usc);

for k=1:length(outTraces)
    if outTraces(k).sampleCount > minSc
        %chop off the rest of the samples (which should be only one, but
        %I'm making it general anyways
        outTraces(k).data = outTraces(k).data(1:minSc);
        outTraces(k).sampleCount = length(outTraces(k).data);
    end
end

%todo: I should add a 'comment'