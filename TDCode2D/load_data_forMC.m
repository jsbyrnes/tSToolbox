function [ dataStruct ] = load_data_forMC( TD_parameters )

    allLats   = [];
    allLons   = [];
    allTS     = [];
    dataE     = [];
    allSta    = {};

    fnames = dir([ TD_parameters.data_dir '/*Measurement.mat']);    
    
    for k = 1:length(fnames)
    
        fname = [ TD_parameters.data_dir '/' fnames(k).name];

        load(fname, 'ts_run', 'Traces');
        traces=ts_run; clear ts_run;

        if TD_parameters.QC

            passedQC = [Traces.QC];
            traces   = traces(passedQC);

        else

            clear Traces

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

        for i = 1:length(TD_parameters.stations_to_remove)

            indkill = strcmp(stations_to_remove(i), sta);

            sta(indkill)      = [];
            lats(indkill)     = [];
            lons(indkill)     = [];
            tS(indkill)       = [];

        end

        allLats  = [allLats lats];
        allLons  = [allLons lons];
        allTS    = [allTS tS];
        allSta   = [allSta sta];
        dataE    = [dataE lats*0+k]; % this is just the event number
    
    end
    
    uSta=unique(allSta); %list of all station names

    allSig  = 0.3*ones(size(allTS));%dummy value

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

