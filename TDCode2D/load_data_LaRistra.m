function [ dataStruct ] = load_data_LaRistra( TD_parameters )

    allLats   = [];
    allLons   = [];
    allTS     = [];
    allSig    = [];
    dataE     = [];
    allSta    = {};

    stations_to_remove       = { 'MB01' 'MB04' 'MB04B' 'MB05' };
    
    fnames = dir('./LaRistraFirstPick/*Measurement.mat');
    fnames = { fnames.name };

    for k=1:length(fnames)

        fname=[ './LaRistraFirstPick/' fnames{k} ];

        load(fname, 'ts_run', 'Traces');

        use = [Traces.QC];
        ts_run = ts_run(use);

        %loop to remove doubled lats/lons (not my preferred solution)
        for kt=1:length(ts_run)

            lt=ts_run(kt).latitude;
            ts_run(kt).latitude=lt(1);

            ln=ts_run(kt).longitude;
            ts_run(kt).longitude=ln(1);

        end

        sta   = {ts_run.station};
        lats  = [ts_run.latitude];
        lons  = [ts_run.longitude];
        tS    = [ts_run.tStar_WF];
        sig   = 0.3*ones(size([ts_run.tStar_WF]));%[traces.sigma];

        for i = 1:length(stations_to_remove)

            indkill = strcmp(stations_to_remove(i), sta);

            sta(indkill)      = [];
            lats(indkill)     = [];
            lons(indkill)     = [];
            tS(indkill)       = [];
            %error(indkill)    = [];

        end
        
        allLats  = [ allLats lats ];
        allLons  = [ allLons lons ];
        allTS    = [ allTS tS ];
        allSig   = [ allSig sig ];
        allSta   = [ allSta sta ];
        dataE    = [ dataE lats*0+k ]; % this is just the event number

    end

    uSta=unique(allSta); %list of all station names

    %%%%synthetic data
    %synthetic model with two steps
    % allTS(:)                                   = 0;%normrnd(0, noise, size(allTS));
    % allTS( allLons < -79 & allLons > -80 )     =  0.2;
    % allTS( allLons < -81 & allLons > -82 )     = -0.2;
    % 
    % for k = 1:3
    % 
    % 	allTS(dataE == k) = allTS(dataE==k) + normrnd(0,0.075, size(allTS(dataE==k)));
    % 
    % end
    % 
    % for k = 4:6
    % 
    % 	allTS(dataE == k) = allTS(dataE==k) + normrnd(0, 0.15, size(allTS(dataE==k)));
    % 
    % end

    % 
    % trueModel = zeros(size(allLons));
    % trueModel( allLons < -81 & allLons > -82 )  = -0.2;
    % trueModel( allLons < -79 & allLons > -80 )      = 0.2;

    %synthetic model with a linear increase
    %allLons = linspace(min(allLons), max(allLons), length(allLons));
    %allLats = linspace(min(allLats), max(allLats), length(allLats));
    %allTS(:) = normrnd(0, noise, size(allTS)) + (( 0.4/(max(allLons) - min(allLons)) *(allLons - min(allLons))) - 0.2);

    %make half the data random
    %allTS( (length(allTS)+1):2*length(allTS) ) = 0.8*rand(size(allTS)) - 0.4;
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
    %dataY = zeros(size(dataY));

    %make the grid on which to invert

    minX = min(dataX) - TD_parameters.buffer;
    maxX = max(dataX) + TD_parameters.buffer;
    minY = min(dataY) - TD_parameters.buffer;
    maxY = max(dataY) + TD_parameters.buffer;

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
    
end

