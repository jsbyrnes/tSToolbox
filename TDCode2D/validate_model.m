function [model] = validate_model(model, dataStruct, TD_parameters)

    for k = 1:length(dataStruct.dataX)
       
        [ ~, closest(k) ] = min( (dataStruct.dataX(k) - model.xCell).^2 + (dataStruct.dataY(k) - model.yCell).^2 );
        
    end

    for k = 1:length(model.xCell)
       
        total(k) = sum(closest == k);
        
    end

    if any(total < TD_parameters.min_sta_per_cell)
        
        model.valid = 0;
        
    else
        
        model.valid = 1;
        
    end
    
end

