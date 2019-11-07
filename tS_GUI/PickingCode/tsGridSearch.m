function [bestAttWf] = tsGridSearch(obsWf,synthWf,tStarV)
% here, pass an observed waveform, a synthetic waveform and a range of
% tstars
% obsWf and synthWf shoud be in irisFetch 'trace' format?
% the synthetic waveform is the source estimate
% out comes the best-fit attenuated source as a structure. In that
% structure there are fields for the tStar and the misfit.

%tStar = the best fit t*
%tStarV = the vector containingthe values f t* t try
%tS = the t* value being tried


attWfs=repmat(synthWf,[length(tStarV),1]); % array of structure to hold the attenuated waveforms

[dummy, btmark]=max(abs(obsWf.data)); %last sample to consider for alignment (big through pick)
% There are alternatives to the last line
% that last line may cause problems.


tic


for k=1:length(tStarV)

    tS=tStarV(k);
   
    
    %t* dependent width of E3: (see widths_calc.m)
    %this is used to calculate lastsamp and in the weighting
    wth=round( min([20 (10+tS*10)]) ); % JUST GUESSED AT SOME NUMBERS
    hwth=round(wth/2); %hte half width
    lastsamp=btmark+hwth;
    %in case the trace is weird and lastsamp is right by the end
    lastsamp=min(lastsamp,length(obsWf.data));
    %it's terrible that this is all hardwired
    
    %attenuate and re-normalize
    synthWfA=wfAttenuate(synthWf,tS); 
    synthWfA.data=synthWfA.data/max(abs(synthWfA.data)); %Attenuated source
    %A for attenuated.
    
    %align using waveform difference
    synthWfAA=synthWfA;
    synthWfAA.data=linemup(synthWfA.data,obsWf.data,lastsamp); %<<<<<<<<<< here probably replace
    %AA for attenuated and Aligned
    
    
    %weighting the misfit
    d=synthWfAA.data(100:lastsamp)-obsWf.data(100:lastsamp); %this is the window
    
%     %this is for checking the time window used for calculating misfit
%     figure(56); clf; hold on
%     plot(synthWfAA.data,'k-')
%     plot(obsWf.data)
%     plot([100,lastsamp],[0 0],'m*')
%     pause
    
    wghts=ones(size(d));
%     wghts(50:ltmark-400)=3; %initial part
%     wghts( ltmark-400 : (btmark-400)-hwth )=.1; 
%     wghts=wghts/3;
    
    d=d.*wghts;
    
    attWfs(k).data=synthWfAA.data; 
    %attWfs stores the best-aligned attenuated source.
    %Later the very best-fitting one is chosen from these
    attWfs(k).misfit=norm(d);
    attWfs(k).tStar=tS;
  
    
end

misfits=[attWfs.misfit];

[berr, bk]=min(misfits); %find which tS produced the smallest misfit

bestAttWf=attWfs(bk);






























