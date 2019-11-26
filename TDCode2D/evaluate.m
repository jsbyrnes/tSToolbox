function [ phiM, likelyhood ] = evaluate( xCell, yCell, tScell, allTS, allSig, dataE, dataX, dataY, TD_parameters)

    F    = scatteredInterpolant(xCell(:), yCell(:), tScell(:), TD_parameters.interp_style, TD_parameters.interp_style);
    ptS = F(dataX,dataY);

    if isempty(ptS)
        
        phiM       = nan;
        likelyhood = nan;
        return
        
    end
    
    if any(any(isnan(ptS)))%nans show up for very closely spaced nodes (maybe?)

        F    = scatteredInterpolant(xCell(:), yCell(:), tScell(:), TD_parameters.interp_style, TD_parameters.interp_style);
        ptS(isnan(ptS)) = F(dataX(isnan(ptS)),dataY(isnan(ptS)));

    end

    %re-demean. The model may have nonzero mean but that's fine
    for k = 1:length(unique(dataE))

        ptS(dataE==k) = ptS(dataE==k) - mean(ptS(dataE==k));

    end

    if strcmp(TD_parameters.likelyhood, 'Gaussian')
    
        phiM       = ((ptS - allTS))*diag(((1./allSig).^2))*((ptS - allTS))';
        likelyhood = -sum(log(allSig*sqrt(2*pi)))*length(allTS) - 0.5*sum(((ptS - allTS)./allSig).^2);
        
    elseif strcmp(TD_parameters.likelyhood, 'Laplacian')
        
        phiM       = sum(abs(ptS - allTS)./allSig);
        likelyhood = -sum(log(2*allSig))*length(allTS) - sum(phiM);
        
    end
            
end