%fit a smooth surface to all t* observations through an inversion scheme.
clear, close all

TD_parameters = define_TDstructure( );
dataStruct    = load_data_forMC(TD_parameters);

%load('BayesianInv_HighLavaPlain.mat', 'dataStruct');

%% Let's invert
parfor k = 1:TD_parameters.n_chains
    
    disp(['Chain #' num2str(k) ]);
    warning off MATLAB:scatteredInterpolant:DupPtsAvValuesWarnId
    models(:, k) = TD_inversion_function(TD_parameters, dataStruct, k);
        
end

disp(['Saving models named ' TD_parameters.savename]);
save(TD_parameters.savename)

models          = models(:);
[xMat, yMat]    = meshgrid(dataStruct.xVec, dataStruct.yVec);
[ model_stats ] = make_map_models( models, xMat, yMat, TD_parameters.sig_flag, -0.5:0.05:0.5, 0.66, TD_parameters.interp_style);

figure
contourf(xMat, yMat, model_stats.mean, -0.15:0.02:0.15);
h = colorbar; h.Label.String = '\Deltat*';
colormap('jet')
hold on
plot(dataStruct.dataX, dataStruct.dataY, 'rs');
