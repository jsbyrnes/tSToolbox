function [ dataStruct ] = load_data_SaltonTrough( TD_parameters )

    allLats   = [];
    allLons   = [];
    allTS     = [];
    allSig    = [];
    dataE     = [];
    allSta    = {};

    %fnames = dir('../DeepResults/*result.mat');
    %fnames = dir('./*witherror.mat');
    %fnames = dir('./*_ShallowWithError.mat');
    %fnames = dir('./tSMeasurements/tS-*');

    dir                      = 'C:\Research\tstar\SaltonTrough\SaltonTrough\CalTech_firstpass\'; 
    dir_ssip                 = 'C:\Research\tstar\SaltonTrough\SaltonTrough\Post2pifix_measurements\';
    QC                       = 1;
    fnames                   = { '004', '006', '011', '020', '028', '030', '046', '048' };
    eventind                 = [ 4 6 11 20 28 30 46 48 ];
    stations_to_remove       = { 'MLAC' 'SNCC' 'SCI2' 'PDM' 'B01' 'M02' 'S01' 'N01'};
    ssip_aswell              = 1;
    pert_ssip                = 0.0;

    load('C:\Research\tstar\SaltonTrough\SaltonTrough\SSIP_data\SSIP.mat');
        
    for k=1:length(fnames)

%         fname=[ dir 'CalTech_SSIP_' fnames{k} '_3Hz_result.mat'];
%         %disp(fname)
% 
%         load(fname, 'ts_run');
%         traces=ts_run; clear ts_run;
% 
%         if QC
% 
%             fname    = [dir 'CalTech_SSIP_' fnames{k} '_3Hz_resultQC.mat'];
%             load(fname);
%             passedQC = [Traces.QC];
%             traces   = traces(passedQC);
%             %error    = error(passedQC);
% 
%         end
% 
%         %loop to remove doubled lats/lons (not my preferred solution)
%         for kt=1:length(traces)
% 
%             lt=traces(kt).latitude;
%             traces(kt).latitude=lt(1);
% 
%             ln=traces(kt).longitude;
%             traces(kt).longitude=ln(1);
% 
%         end
% 
%         sta   = {traces.station};
%         lats  = [traces.latitude];
%         lons  = [traces.longitude];
%         tS    = [traces.tStar];
% 
%         for i = 1:length(stations_to_remove)
% 
%             indkill = strcmp(stations_to_remove(i), sta);
% 
%             sta(indkill)      = [];
%             lats(indkill)     = [];
%             lons(indkill)     = [];
%             tS(indkill)       = [];
%             %error(indkill)    = [];
% 
%         end

        if ssip_aswell

            fname=[ dir_ssip 'SSIP_' fnames{k} '_1.5Hz_result.mat'];

            load(fname, 'ts_run');
            traces=ts_run; clear ts_run;

            if QC

                fname=[ dir_ssip 'SSIP_' fnames{k} '_1.5Hz_resultQC.mat'];
                load(fname);
                passedQC = [Traces.QC];
                traces   = traces(passedQC);
                %error    = error(passedQC);

            end

            %loop to remove doubled lats/lons (not my preferred solution)
            for kt=1:length(traces)

                lt=traces(kt).latitude;
                traces(kt).latitude=lt(1);

                ln=traces(kt).longitude;
                traces(kt).longitude=ln(1);

            end

            sta   = {traces.station};
            lats  = [traces.latitude];
            lons  = [traces.longitude];
            tS    = [traces.tStar];

            for i = 1:length(stations_to_remove)

                indkill = strcmp(stations_to_remove(i), sta);

                sta(indkill)      = [];
                lats(indkill)     = [];
                lons(indkill)     = [];
                tS(indkill)       = [];
                %error(indkill)    = [];

            end

%             sta   = [ sta sta ];
%             lats  = [ lats lats ];
%             lons  = [ lons lons ];
%             tS    = [ tS (tS + pert_ssip) ];
    %         sta   = [ sta2 ];
    %         lats  = [ lats2 ];
    %         lons  = [ lons2 ];
    %         tS    = [ (tS2 + pert_ssip) ];

        end

        allLats  = [allLats lats];
        allLons  = [allLons lons];
        allTS    = [allTS tS];
        allSta   = [allSta sta];
        dataE    = [dataE lats*0+k]; % this is just the event number
        %allError = [allError error]
        
    end

    uSta=unique(allSta); %list of all station names

    %%%%synthetic data
    % allTS(:)                                   = normrnd(0, noise, size(allTS));
    % allTS( allLons < -79 & allLons > -80 )     = 0.2 + normrnd(0, noise, size(allTS( allLons < -79 & allLons > -80 )));
    % allTS( allLons < -81 & allLons > -82 )     = -0.2 + normrnd(0, noise, size(allTS( allLons < -81 & allLons > -82 )));
    % allTS( (length(allTS)+1):2*length(allTS) ) = 0.8*rand(size(allTS)) - 0.4;
    allSig  = 0.3*ones(size(allTS));
    % dataE   = [ dataE ( dataE + max(dataE)) ];
    % allLats = [ allLats allLats ];
    % allLons = [ allLons allLons ];
    % allSta  = [ allSta allSta ];

    for k = 1:length(unique(dataE))

        allTS(dataE==k) = allTS(dataE==k) - mean(allTS(dataE==k));

    end

    %% define a map projection and inversion grid
    %origin in the center of the array footprint

    centerLat = min(allLats)+(max(allLats)-min(allLats))/2;
    centerLon = min(allLons)+(max(allLons)-min(allLons))/2;

    origin               = [ centerLat centerLon ];           % array center [ lat lon ]
    mstruct              = defaultm('mercator');
    mstruct.origin       = [ origin TD_parameters.rotation ];%second number is rotation
    mstruct              = defaultm( mstruct );
    mstruct.scalefactor  = 6371;

    %%%%%%%%%%
    %%%consider rotating the grid for the inversion so cut down the number
    %%%of unsed nodes
    %%%%%%%%%%

    %projected lat/lon of observations to X/Y:
    [dataX, dataY] = mfwdtran(mstruct,allLats,allLons);
    %figure; plot(dataX,dataY,'k.');
    %dataY(:) = 0;
    %make the grid on which to invert

    load coastlines.mat

    [coastX, coastY] = mfwdtran(mstruct, coastlat, coastlon);
    
    minX = min(dataX) - TD_parameters.buffr;
    maxX = max(dataX) + TD_parameters.buffr;
    minY = min(dataY) - TD_parameters.buffr;
    maxY = max(dataY) + TD_parameters.buffr;

    xVec = minX:TD_parameters.nodeSpacing:maxX;
    yVec = minY:TD_parameters.nodeSpacing:maxY;    
    
    dataStruct.allTS   = allTS;
    dataStruct.allLats = allLats;
    dataStruct.allLons = allLons;
    dataStruct.allSig  = allSig;
    dataStruct.allSta  = allSta;
    dataStruct.dataX   = dataX;
    dataStruct.dataY   = dataY;
    dataStruct.dataE   = dataE;
    dataStruct.xVec    = xVec;
    dataStruct.yVec    = yVec;
    dataStruct.coastX  = coastX;
    dataStruct.coastY  = coastY;
            
end

