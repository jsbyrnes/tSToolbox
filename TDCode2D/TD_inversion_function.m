function [ model_hist ] = TD_inversion_function(TD_parameters, dataStruct)

    rng(round(mod(now*1e12,1e3)))%seed down to less than the milisecond

    xVec = dataStruct.xVec;
    yVec = dataStruct.yVec;

    n = length(dataStruct.dataX);
               
    %initialize the chain
    %first four nodes are top corners and CANNOT be moved
    
    model.valid = 0;
    
    while ~model.valid
    
        model.nCells          = randi([ TD_parameters.min_cells TD_parameters.max_cells ]);
        model.xCell           = min(xVec(:)) + (max(xVec(:)) - min(xVec))*rand([ model.nCells 1 ]);
        model.yCell           = min(yVec(:)) + (max(yVec(:)) - min(yVec))*rand([ model.nCells 1 ]);
    
        model = validate_model(model, dataStruct, TD_parameters);
        
    end
    
    model.tSCell  = (TD_parameters.range(2) - TD_parameters.range(1))*rand([ model.nCells 1 ]) + TD_parameters.range(1);

    %if you want to see it, plot with the voronoi command
    %or interpolate to the grid
    
    %initalize to the mean sigma for the event, and set a dSig for an events
    %deviation from that mean. Then perturb the mean sigma only. Only used if
    %sig_flag is set to 2
    dX          = unique(dataStruct.dataX);

    if TD_parameters.sig_flag == 1
        
        allSig = dataStruct.allSig;    
        %model.allSig = rand(1)*allSig(1);
        model.allSig = allSig(1);
        allSig = model.allSig*ones(size(allSig));    
        
    elseif TD_parameters.sig_flag == 2

        station_sig = zeros(size(dX));
        for k = 1:length(dX)

            station_sig(k) = rand(1)*mean(dataStruct.allSig( dataStruct.dataX == dX(k)));

        end
        model.allSig = station_sig;

        for k = 1:length(dataStruct.allTS)
        
            allSig(k) = station_sig(dataStruct.dataX(k) == dX(k));
            
        end
        
    elseif TD_parameters.sig_flag == 3

        evt_sig     = zeros(max(dataStruct.dataE),1);
        for k = 1:max(dataStruct.dataE)

            evt_sig(k) = rand(1)*mean(dataStruct.allSig( dataStruct.dataE == k ));

        end
        model.allSig = evt_sig;
        
        for k = 1:length(dataStruct.allTS)
        
            allSig(k) = evt_sig(dataStruct.dataE(k));
            
        end
                
    else
        
        error('fix sig_flag');

    end
    
    model.phi = NaN;
        
    [ model.phi, likelyhood] = evaluate(model.xCell, model.yCell, model.tSCell, dataStruct.allTS, allSig, ...
        dataStruct.dataE, dataStruct.dataX, dataStruct.dataY, TD_parameters);
            
    model.likelyhood         = likelyhood;
        
    for iter = 1:TD_parameters.n_iter
        
        %pick birth, death, move, change, or noise
        if TD_parameters.sig_sig > 0

            action = randi([1 5],1);

        else

            action = randi([1 4],1);

        end

        model.action = action;
        model.accept = 0;

        switch action

            case 1 %birth
                
                valid = 0;
                
                while ~valid
                
                    xNew   = rand(1)*(max(xVec) - min(xVec)) + min(xVec);%xVec(randi(length(xVec),[1 1]));
                    yNew   = rand(1)*(max(yVec) - min(yVec)) + min(yVec);

                    %get the current tS for this point and make a new one
                    F     = scatteredInterpolant(model.xCell, model.yCell, model.tSCell, TD_parameters.interp_style, 'linear');
                    ctS   = F(xNew, yNew);
                    tSnew = normrnd(ctS, TD_parameters.sig_b);

                    xCell_n  = [ model.xCell;  xNew ];
                    yCell_n  = [ model.yCell;  yNew ];
                    tSCell_n = [ model.tSCell; tSnew ];

                    modelN      = model;
                    modelN.xCell = xCell_n;
                    modelN.yCell = yCell_n;
                    
                    modelN = validate_model(modelN, dataStruct, TD_parameters);
                    
                    valid = modelN.valid;
                    
                end
                
                tSCell_n(end) = max( [ TD_parameters.range(1) min([ tSCell_n(end) TD_parameters.range(2) ]) ]);
                
                [phiN, likelyhood] = evaluate(xCell_n, yCell_n, tSCell_n, dataStruct.allTS, allSig, dataStruct.dataE, ...
                    dataStruct.dataX, dataStruct.dataY, TD_parameters);
                
                if strcmp(TD_parameters.likelyhood, 'Gaussian')
                    
                    alpha = min([1 (TD_parameters.sig_b*sqrt(2*pi))/diff(TD_parameters.range)*...
                        exp( ((tSnew - ctS)^2)/(2*TD_parameters.sig_b^2) - (phiN - model.phi)/2)]);
                    
                elseif strcmp(TD_parameters.likelyhood, 'Laplacian')
                    
                    alpha = min([1 (TD_parameters.sig_b*sqrt(2*pi))/diff(TD_parameters.range)*...
                        exp( ((tSnew - ctS)^2)/(2*TD_parameters.sig_b^2) - (phiN - model.phi))]);
                    
                end
                
                if rand(1) <= alpha && model.nCells < TD_parameters.max_cells
                    
                    model.nCells        = model.nCells + 1;
                    model.xCell         = xCell_n;
                    model.yCell         = yCell_n;
                    model.tSCell        = tSCell_n;
                    model.phi           = phiN;
                    model.likelyhood    = likelyhood;
                    model.accept        = action;
                    
                end

            case 2 %death. No need to validate

                if model.nCells > TD_parameters.min_cells

                    kill = randi([ 5 model.nCells ], 1);

                    xCell_n        = model.xCell;
                    yCell_n        = model.yCell;
                    tSCell_n       = model.tSCell;
                    ctS            = tSCell_n(kill);
                    
                    tSCell_n(kill) = [];
                    xCell_n(kill)  = [];
                    yCell_n(kill)  = [];

                    F       = scatteredInterpolant(xCell_n, yCell_n, tSCell_n, TD_parameters.interp_style, 'linear');
                    newtS   = F(model.xCell(kill), model.yCell(kill)); 
                    
                    [phiN, likelyhood] = evaluate(xCell_n, yCell_n, tSCell_n, dataStruct.allTS, allSig, ...
                        dataStruct.dataE, dataStruct.dataX, dataStruct.dataY, TD_parameters);
                                        
                    if strcmp(TD_parameters.likelyhood, 'Gaussian')

                        alpha = min([1 diff(TD_parameters.range)/(TD_parameters.sig_b*sqrt(2*pi))*...
                            exp( -((ctS - newtS)^2)/(2*TD_parameters.sig_b^2) - (phiN - model.phi)/2)]);

                    elseif strcmp(TD_parameters.likelyhood, 'Laplacian')

                        alpha = min([1 diff(TD_parameters.range)/(TD_parameters.sig_b*sqrt(2*pi))*...
                            exp( -((ctS - newtS)^2)/(2*TD_parameters.sig_b^2) - (phiN - model.phi))]);

                    end
                        
                    if rand(1) <= alpha

                        model.nCells        = model.nCells - 1;
                        model.xCell         = xCell_n;
                        model.yCell         = yCell_n;
                        model.tSCell        = tSCell_n;
                        model.phi           = phiN;
                        model.likelyhood    = likelyhood;
                        model.accept        = action;

                    end

                end

            case 3 %change t* value. No need to validate. 
                
                change           = randi(model.nCells, 1);
                tSCell_n         = model.tSCell;
                tSCell_n(change) = normrnd(model.tSCell(change), TD_parameters.sig_b);

                tSCell_n(change) = max( [ TD_parameters.range(1) min([ tSCell_n(change) TD_parameters.range(2) ]) ]);

                [phiN, likelyhood] = evaluate(model.xCell, model.yCell, tSCell_n, dataStruct.allTS, allSig, dataStruct.dataE, ...
                    dataStruct.dataX, dataStruct.dataY, TD_parameters);            
                
                r = rand(1);

                if strcmp(TD_parameters.likelyhood, 'Gaussian')
                    
                    alpha = min( [ 1 exp(-(phiN - model.phi)/2) ] );
                    
                elseif strcmp(TD_parameters.likelyhood, 'Laplacian')
                    
                    alpha = min( [ 1 exp(-(phiN - model.phi)) ] );
                    
                end
                
                if r <= alpha

                    model.tSCell        = tSCell_n;
                    model.phi           = phiN;
                    model.likelyhood    = likelyhood;
                    model.accept        = action;

                end

            case 4 %move

                valid = 0;
                
                while ~valid
                
                    move    = randi(model.nCells, 1);
                    xCell_n = model.xCell;
                    yCell_n = model.yCell;

                    xCell_n(move) = normrnd(model.xCell(move), (max(xVec) - min(xVec))*(TD_parameters.sig_r/100));
                    yCell_n(move) = normrnd(model.yCell(move), (max(yVec) - min(yVec))*(TD_parameters.sig_r/100))/TD_parameters.ydamp;
                    
                    if xCell_n(move) < max(xVec) && xCell_n(move) > min(xVec) ...
                            && yCell_n(move) < max(yVec) && yCell_n(move) > min (yVec)
                    
                        valid = 1;
                        
                    end
                    
                    if valid
                    
                        modelN = model;
                        modelN.xCell = xCell_n;
                        modelN.yCell = yCell_n;

                        modelN = validate_model(modelN, dataStruct, TD_parameters);

                        valid = modelN.valid;
                    
                    end
                        
                end
                
                [phiN, likelyhood] = evaluate(xCell_n, yCell_n, model.tSCell, dataStruct.allTS, allSig, dataStruct.dataE, ...
                    dataStruct.dataX, dataStruct.dataY, TD_parameters);

                r = rand(1);
                
                if strcmp(TD_parameters.likelyhood, 'Gaussian')
                    
                    alpha = min( [ 1 exp(-(phiN - model.phi)/2) ] );
                    
                elseif strcmp(TD_parameters.likelyhood, 'Laplacian')
                    
                    alpha = min( [ 1 exp(-(phiN - model.phi)) ] );
                    
                end

                if r <= alpha

                    model.xCell         = xCell_n;
                    model.yCell         = yCell_n;
                    model.phi           = phiN;
                    model.likelyhood    = likelyhood;
                    model.accept        = action;

                end

            case 5 %change sigma
                
                if TD_parameters.sig_flag == 1 %uniform error
                    
                    sig_n = -1;
                    while sig_n <= 0
                        
                        sig_n       = normrnd(allSig(1), TD_parameters.sig_sig);%assumes uniform
                        
                    end
                    
                    [phiN, likelyhood] = evaluate(model.xCell, model.yCell, model.tSCell, dataStruct.allTS, sig_n*ones(size(allSig)), ...
                        dataStruct.dataE, dataStruct.dataX, dataStruct.dataY, TD_parameters);
                    
                    r = rand(1);
                    
                    if strcmp(TD_parameters.likelyhood, 'Gaussian')

                        alpha = min( [ log(1) (log((model.allSig/sig_n))*(n) - ((phiN - model.phi)/2)) ] );

                    elseif strcmp(TD_parameters.likelyhood, 'Laplacian')

                        alpha = min( [ log(1) (log((model.allSig/sig_n))*(n) - (phiN - model.phi)) ] );

                    end
                    
                    if log(r) <= alpha
                        
                        allSig              = sig_n*ones(size(allSig));
                        model.allSig        = sig_n(1);
                        model.phi           = phiN;
                        model.likelyhood    = likelyhood;
                        model.accept        = action;
                        
                    end
                    
                elseif TD_parameters.sig_flag == 2 %per station errors
                    
                    change                = randi(length(station_sig), 1);
                    station_sig_n         = station_sig;
                    
                    station_sig_n(change) = -1;
                    while station_sig_n(change) <= 0
                    
                        station_sig_n(change) = normrnd(station_sig(change), TD_parameters.sig_sig);
                    
                    end
                    allSig_n              = allSig;
                    
                    %update station sigmas
                    for k = 1:length(allSig_n)
                        
                        allSig_n(k) = station_sig_n(dataStruct.dataX(k) == dX);
                        
                    end
                    
                    [phiN, likelyhood] = evaluate(model.xCell, model.yCell, model.tSCell, dataStruct.allTS, allSig_n, dataStruct.dataE, ...
                        dataStruct.dataX, dataStruct.dataY, TD_parameters);
                    
                    r = rand(1);                                
                    if strcmp(TD_parameters.likelyhood, 'Gaussian')

                        alpha = min( [ log(1) (sum(log(allSig)) - sum(log(allSig_n)) - ((phiN - model.phi)/2)) ] );

                    elseif strcmp(TD_parameters.likelyhood, 'Laplacian')

                        alpha = min( [ log(1) (sum(log(allSig)) - sum(log(allSig_n)) - ((phiN - model.phi))) ] );

                    end
                    
                    if log(r) <= alpha
                        
                        allSig              = allSig_n;
                        station_sig         = station_sig_n;
                        model.allSig        = station_sig_n;
                        model.phi           = phiN;
                        model.likelyhood    = likelyhood;
                        model.accept        = action;
                        
                    end
                    
                elseif TD_parameters.sig_flag == 3 %per event errors
                    
                    change            = randi(length(evt_sig), 1);
                    evt_sig_n         = evt_sig;
                    
                    evt_sig_n(change) = -1;
                    while evt_sig_n(change) <= 0
                    
                        evt_sig_n(change) = normrnd(evt_sig(change), TD_parameters.sig_sig);
                        
                    end
                        
                    allSig_n          = allSig;
                    
                    %update station sigmas
                    for k = 1:length(allSig_n)
                        
                        allSig_n(k) = evt_sig_n(dataStruct.dataE(k));
                        
                    end
                    
                    [phiN, likelyhood] = evaluate(model.xCell, model.yCell, model.tSCell, dataStruct.allTS, allSig_n, ...
                        dataStruct.dataE, dataStruct.dataX, dataStruct.dataY, TD_parameters);
                    
                    r = rand(1);
                    if strcmp(TD_parameters.likelyhood, 'Gaussian')

                        alpha = min( [ log(1) (sum(log(allSig)) - sum(log(allSig_n)) - ((phiN - model.phi)/2)) ] );

                    elseif strcmp(TD_parameters.likelyhood, 'Laplacian')

                        alpha = min( [ log(1) (sum(log(allSig)) - sum(log(allSig_n)) - ((phiN - model.phi))) ] );

                    end
                    
                    if log(r) <= alpha
                        
                        allSig              = allSig_n;
                        evt_sig             = evt_sig_n;
                        model.allSig        = evt_sig_n;
                        model.phi           = phiN;
                        model.likelyhood    = likelyhood;
                        model.accept        = action;
                        
                    end
                    
                end
                                                        
        end

        llhist(iter)        = likelyhood;
        model.likelyhood    = likelyhood;
        sig_hist(iter, :)   = model.allSig;
        n_hist(iter)        = model.nCells;
                
        if iter == 1
           
            model_hist = model;
            
        end
        
        if iter >= TD_parameters.burn_in && mod(iter, TD_parameters.keep_each) == 0
        
            model_hist(end + 1) = model;
                        
        end
        
    end

    model_hist(1) = [];

end

