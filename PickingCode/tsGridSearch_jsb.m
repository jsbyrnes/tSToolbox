function [bestAttWf, misfits] = tsGridSearch_jsb(obsWf, synthWf, tStarV, fitting_window, filter_bounds)
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

%editted out by jsb
%[dummy, btmark]=max(abs(obsWf.data)); %last sample to consider for alignment (big through pick)
% There are alternatives to the last line
% that last line may cause problems.

owf_data = obsWf.data;

for k=1:length(tStarV)

    tS=tStarV(k);
           
    %attenuate and re-normalize
    synthWfA=wfAttenuate(synthWf,tS, 0.001); 
    synthWfA.data=synthWfA.data/max(abs(synthWfA.data)); %Attenuated source
    %A for attenuated.
    
    %align using waveform difference
    synthWfAA=synthWfA;
    
    %added by JSB if there is a filter, you need to apply it to the
    %synthetic
%     if ~isnan(filter_bounds)
%         
%         synthWfAA = wfButterworth( synthWfAA, filter_bounds);
% 
%     end
    
    %add by jsb to address normalization issues
    %normalization
%     obsWf.data       = obsWf.data/(sum(abs(obsWf.data(fitting_window(1):fitting_window(2)))));
%     synthWfAA.data   = synthWfAA.data/(sum(abs(synthWfAA.data(fitting_window(1):fitting_window(2)))));
%     obsWf.data       = obsWf.data/(sum(abs(obsWf.data)));
%     synthWfAA.data   = synthWfAA.data/(sum(abs(synthWfAA.data)));

    owf_data         = owf_data/rms(owf_data(fitting_window(1):fitting_window(2)));
    synthWfAA.data   = synthWfAA.data/rms(synthWfAA.data(fitting_window(1):fitting_window(2)));
    
    %linemup pads the synthetic time series, and that is throughing the
    %fitting window off. Save the length of the window and adjust
    
    %jsb changed to use the window that the user selected
    
    if length(synthWfAA.data) > length(owf_data)
        
        synthWfAA.data = synthWfAA.data(1:length(owf_data));
        
    elseif length(synthWfAA.data) < length(owf_data)
        
        owf_data = owf_data(1:length(synthWfAA.data));
        
    end
    
    [ synthWfAA.data ] = linemup_quick(synthWfAA.data, owf_data);
    
    d = synthWfAA.data(fitting_window(1):fitting_window(2)) ...
                - owf_data(fitting_window(1):fitting_window(2)); %this is the window
                        
    attWfs(k).data=synthWfAA.data; 
    %attWfs stores the best-aligned attenuated source.
    %Later the very best-fitting one is chosen from these
    attWfs(k).misfit = rms(d);
    attWfs(k).tStar_WF=tS;
    
end

misfits=[attWfs.misfit];

[berr, bk]=min(misfits); %find which tS produced the smallest misfit

% figure
% plot(tStarV, misfits)
% xlabel('t*')
% ylabel('Misfits')

bestAttWf=attWfs(bk);

%bestAttWf.tStar_FD = spectral_slopes_tS(bestAttWf.sampleRate, filter_bounds, obsWf.data, synthWf.data, []);





























