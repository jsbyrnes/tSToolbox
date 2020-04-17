function [ TD_parameters ] = define_TDstructure( )

    %%%%%%%search parameters
    TD_parameters.sig_b            = 0.01;%in s
    TD_parameters.sig_r            = 20;
    TD_parameters.sig_sig          = 0.001;%in s
    TD_parameters.sig_flag         = 3;%1, uniform. 2, station. 3, event.
    TD_parameters.min_error        = 0.0001;
    TD_parameters.interp_style     = 'nearest';%linear or nearest
    TD_parameters.ydamp            = 1;
    TD_parameters.range            = [ -0.1 0.1 ];
    TD_parameters.n_chains         = 48;
    TD_parameters.n_iter           = 5e5;%for now, just iterate to a max
    TD_parameters.burn_in          = 2.5e5;
    TD_parameters.keep_each        = 2.5e3;
    TD_parameters.max_cells        = 100;%for starting
    TD_parameters.min_cells        = 5;%for entire model, CANNOT be below 4
    TD_parameters.min_sta_per_cell = 2; %restricts the minimum scale features in the model. 
                                        %otherwise the code tends to carve
                                        %out weird stations and give them
                                        %their own cells
    TD_parameters.likelyhood       = 'Laplacian';%Gaussian or Laplacian
    
    %%%%%%%map parameters
    TD_parameters.nodeSpacing      = 20;%also used for slices the ray path
    TD_parameters.buffer           = 50;%in km
    TD_parameters.rotation         = 0;%in degrees
    TD_parameters.save_for_gmt     = 1;

end

