function [ al_tr ] = linemup( swf, owf, lastsamp )
% swf=synthetic waveform, owf=observed waveform

madeSwitch=false;
%make rows if colummns
if size(swf,1) == length(swf)
    madeSwitch=true;
    swf=swf';
    owf=owf';
end


swf(1:100)=[]; % this gets rid of the first 100 samples.

maxlag=200;
linerr=nan(1,maxlag);
atrM=nan(maxlag,length(owf));

for lag=1:maxlag
    pad1=zeros(1,lag);
    pad2=zeros(1,length(owf)-(length(swf)+lag));
    a1=[pad1 swf pad2]; %makes the array at least as long as auxwf
    a1=a1(1:length(owf)); %shortens it if necessary
     
    atrM(lag,:)=a1;
    
    try
    linerr(lag)=norm(a1(100:lastsamp)-owf(100:lastsamp));
    catch
        keyboard
    end
end

[m, shft]=min(linerr);
al_tr=atrM(shft,:);

%you can probably do this without a for loop 


if madeSwitch;
   al_tr=al_tr'; 
end

end

