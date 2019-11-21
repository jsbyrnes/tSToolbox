function [ model_stats ] = make_map_models( models, xMat, yMat, sig_flag, pdf_x, unc_level, interp_style)

    warning off MATLAB:scatteredInterpolant:DupPtsAvValuesWarnId
    
    for k = 1:length(models)
        
        waitbar(k/length(models))
        
        fine_models.model = models(k);
        fine_models.xMat  = xMat;
        fine_models.yMat  = yMat;
        
        F    = scatteredInterpolant(models(k).xCell(:), ...
            models(k).yCell(:), models(k).tSCell(:), interp_style, interp_style);
        fine_models.tS = F(xMat,yMat); 
        
        if any(any(isnan(fine_models.tS)))
            
            F    = scatteredInterpolant(fine_models.xMat(~isnan(fine_models.tS)), ...
                fine_models.yMat(~isnan(fine_models.tS)), fine_models.tS(~isnan(fine_models.tS)), 'nearest', 'nearest');
            fine_models.tS(isnan(fine_models.tS)) = F(xMat(isnan(fine_models.tS)),yMat((isnan(fine_models.tS))));
            
        end
        
        fine_models.tS    = fine_models.tS - mean(fine_models.tS(:));
        
        %how big is each cell? Gives an idea of the model smoothness
        for j = 1:numel(xMat)
           
            cellDistancevec(j) = min( sqrt((xMat(j) - models(k).xCell).^2 + ...
                (yMat(j) - models(k).yCell).^2 ));
                        
        end
                    
        modelset(:,:,k) = fine_models.tS;
                                 
        if sig_flag > 1
            
            sigmat(k, :) = models(k).allSig;
            
        end
        
        clear fine_models cellDistancevec
        
    end
        
    if sig_flag == 1
        
        model_stats.sig     = mean([models.allSig]);
        
    else
            
         model_stats.sig     = mean(sigmat);
         model_stats.sig_sig = std([models.allSig], 0, 2);
         
    end
    
    for k = 1:numel(xMat)
               
        waitbar(k/numel(xMat))
        
        [i,j] = ind2sub(size(xMat), k);
        
        model_stats.mean(k)      = mean(modelset(i,j,:));
        model_stats.std(k)       = std(modelset(i,j,:));
        model_stats.skewness(k)  = skewness(modelset(i,j,:));
        model_stats.kurtosis(k)  = kurtosis(modelset(i,j,:));
        
        pdf_x_tmp = linspace(min(pdf_x), max(pdf_x), length(pdf_x) + 1);
        pdf_set = histcounts(squeeze(modelset(i,j,:)), pdf_x_tmp, 'Normalization', 'pdf');
        %pdf_set = ksdensity(squeeze(modelset(i,j,:)), pdf_x);
        
        pdf_x_fine = linspace(min(pdf_x), max(pdf_x), 10000);
        y = interp1(pdf_x, pdf_set, pdf_x_fine, 'spline');
        
        [~, ind] = max(y);
        model_stats.mode(k) = pdf_x_fine(ind);

        cdf_set(k,:) = cumtrapz(pdf_x, pdf_set);

        [~,ind] = min(abs(cdf_set(k,:) - ( 0.5 - unc_level/2 ) ));
        model_stats.unc_p(k) = pdf_x(ind);
        [~,ind] = min(abs(cdf_set(k,:) - ( 0.5 + unc_level/2 ) ));
        model_stats.unc_n(k) = pdf_x(ind);
            
        pdf_set_plot(k,:) = pdf_set;
        
    end
    
    model_stats.pdf_set   = reshape(pdf_set_plot, [size(xMat) length(pdf_x)]);
    model_stats.mean      = reshape(model_stats.mean, size(xMat));
    model_stats.std       = reshape(model_stats.std, size(xMat));
    model_stats.skewness  = reshape(model_stats.skewness, size(xMat));
    model_stats.kurtosis  = reshape(model_stats.kurtosis, size(xMat));
    model_stats.mode      = reshape(model_stats.mode, size(xMat));
    model_stats.unc_p     = reshape(model_stats.unc_p, size(xMat));
    model_stats.unc_n     = reshape(model_stats.unc_n, size(xMat));
    
end

