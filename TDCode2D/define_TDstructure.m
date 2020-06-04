function [ TD_parameters ] = define_TDstructure( )

    %%%%%%%where to get the data
    TD_parameters.savename           = 'test';
    TD_parameters.data_dir           = 'C:\Research\tstar\RioGrandRift\Rio_nofilt';%path to a folder with *Measurement.mat files
    TD_parameters.QC                 = 1;%remove QCed traces (1) or keep them in (0)
    TD_parameters.stations_to_remove = { };%empty to use all stations
    
    %%%%%%%search parameters
    TD_parameters.sig_b            = 0.01;%in s
    TD_parameters.sig_r            = 20;
    TD_parameters.sig_sig          = 0.001;%in s
    TD_parameters.sig_flag         = 1;%1, uniform. 2, station. 3, event.
    TD_parameters.min_error        = 0.0001;
    TD_parameters.interp_style     = 'nearest';%linear or nearest
    TD_parameters.ydamp            = 0;
    TD_parameters.range            = [ -0.4 0.4 ];
    TD_parameters.n_chains         = 8;
    TD_parameters.n_iter           = 5e5;%for now, just iterate to a max
    TD_parameters.burn_in          = 2.e5;
    TD_parameters.keep_each        = 1e4;
    TD_parameters.print_on         = 0.5;%in percent completed
    TD_parameters.max_cells        = 100;%for starting
    TD_parameters.min_cells        = 5;%for entire model, CANNOT be below 4
    TD_parameters.min_sta_per_cell = 0; %restricts the minimum scale features in the model. Not recommended.

    TD_parameters.likelyhood       = 'Gaussian';%Gaussian or Laplacian
    
    %%%%%%%map parameters
    TD_parameters.nodeSpacing      = 20;
    TD_parameters.buffer           = 50;%in km
    TD_parameters.rotation         = 0;%in degrees

    if TD_parameters.n_iter <= TD_parameters.burn_in
       
        error('Burn-in is larger than the total number of iterations');
        
    end
    
end

