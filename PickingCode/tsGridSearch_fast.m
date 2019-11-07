function [bestAttWf, misfits] = tsGridSearch_fast(obsWf, synthWf, tStarV, btmark, fitting_window, filter_bounds)
% here, pass an observed waveform, a synthetic waveform and a range of
% tstars
% obsWf and synthWf shoud be in irisFetch 'trace' format?
% the synthetic waveform is the source estimate
% out comes the best-fit attenuated source as a structure. In that
% structure there are fields for the tStar and the misfit.

%tStar = the best fit t*
%tStarV = the vector containingthe values f t* t try
%tS = the t* value being tried

%tukeywin on the data
% obsWf.data   = obsWf.data.*tukeywin(length(obsWf.data));
% synthWf.data = synthWf.data.*tukeywin(length(synthWf.data));

attWfs=repmat(synthWf,[length(tStarV),1]); % array of structure to hold the attenuated waveforms

owf_data0 = obsWf.data;

for k=1:length(tStarV)

    tS=tStarV(k);
   
    owf_data = owf_data0;%to get aorund parfor loop
    
    %disp([ 't* of ' num2str(tS) ])
    
    %t* dependent width of E3: (see widths_calc.m)
    %this is used to calculate lastsamp and in the weighting
    wth=round( min([20 (10+tS*10)]) ); % JUST GUESSED AT SOME NUMBERS
    hwth=round(wth/2); %hte half width
    lastsamp=btmark+hwth;
    %in case the trace is weird and lastsamp is right by the end
    lastsamp=min(lastsamp,length(owf_data));
    %it's terrible that this is all hardwired
    
    %attenuate and re-normalize
    synthWfA=wfAttenuate(synthWf,tS, 0.1); 
    synthWfA.data=synthWfA.data/max((synthWfA.data)); %Attenuated source
    %A for attenuated.
    
    %align using waveform difference
    synthWfAA=synthWfA;
    
    %added by JSB if there is a filter, you need to apply it to the
    %synthetic
    if ~isnan(filter_bounds)
        
        synthWfAA = wfCosFilter( synthWfAA, filter_bounds);

    end
    
    %add by jsb to address normalization issues
    %normalization
%     obsWf.data       = obsWf.data/(sum(abs(obsWf.data(fitting_window(1):fitting_window(2)))));
%     synthWfAA.data   = synthWfAA.data/(sum(abs(synthWfAA.data(fitting_window(1):fitting_window(2)))));
%     obsWf.data       = obsWf.data/(sum(abs(obsWf.data)));
%     synthWfAA.data   = synthWfAA.data/(sum(abs(synthWfAA.data)));
    owf_data         = owf_data/max((owf_data));
    synthWfAA.data   = synthWfAA.data/max((synthWfAA.data));
    
    %linemup pads the synthetic time series, and that is throughing the
    %fitting window off. Save the length of the window and adjust
    
    %jsb changed to use the window that the user selected
    synthWfAA.data  = linemup_fast(synthWfAA.data, owf_data); %<<<<<<<<<< here probably replace
    %synthWfAA.data = linemup(synthWfA.data,obsWf.data,lastsamp, fitting_window); %<<<<<<<<<< here probably replace
    %AA for attenuated and Aligned
            
    %renormalize because the fitting window changes
    %JSB update I added this in but I think it is now unnecessary
%     obsWf.data       = obsWf.data/(sum(abs(obsWf.data(window_shift + fitting_window(1):window_shift + fitting_window(2)))));
%     synthWfAA.data   = synthWfAA.data/(sum(abs(synthWfAA.data(window_shift + fitting_window(1):window_shift + fitting_window(2)))));
    
    %weighting the misfit
    
    try
    
        d = synthWfAA.data(fitting_window(1):fitting_window(2)) ...
            - owf_data(fitting_window(1):fitting_window(2)); %this is the window
    
        %JSB for debugging
%         figure(10),clf, plot(obsWf.data), hold on, plot(synthWfAA.data), plot(obsWf.data - synthWfAA.data); ylim([-1 1]);
%         title(num2str(tS))
%         keyboard
        
    catch
        
        keyboard
        
    end
        
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
    attWfs(k).misfit=norm(d)/length(d);
    attWfs(k).tStar=tS;
  
    %jsb debugging code on Oct 11th, 2016
%     figure(2)
%     hold on
%     plot(synthWfAA.data)%, hold on, plot(obsWf.data), plot(synthWfAA.data - obsWf.data(1:length(synthWfAA.data)), 'k');
%     title([ 't* ' num2str(tS) ', misfit ' num2str(attWfs(k).misfit)]);
%     xlim([390 440])
%     %ylim([-4e-3 4e-3])
%     keyboard
    
end

misfits=[attWfs.misfit];

[berr, bk]=min(misfits); %find which tS produced the smallest misfit

bestAttWf=attWfs(bk);





























