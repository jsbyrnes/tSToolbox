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
                    
    elseif strcmp(TD_parameters.likelyhood, 'DoubleGaussian')
        
        phiM       = sum((ptS - allTS).^2);%not used
        likelyhood = sum( model.allSig(3)/(sqrt(2*pi*model.allSig(1)))*...
            exp( ((ptS - allTS).^2)./(2*model.allSig(1)^2)) + ...
            (1 - model.allSig(3))/(sqrt(2*pi*model.allSig(2)))*...
            exp( ((ptS - allTS).^2)./(2*model.allSig(2)^2)));
        
    end
    
end