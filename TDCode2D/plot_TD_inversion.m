%fit a smooth surface to all t* observations through an inversion scheme.
clear
clc
close all
load('CMfine.mat')

savename = 'Rio16';

loadname = './AllRio_tSWF*';

label = '\Deltat*_p, s';
%label = '\Delta\Sigma';

states_or_boundaries = 2;%1 for states, 2 for boundaries
masking_threshold    = 0.1;%make Inf for no masking
cmap                 = flipud(cm);
map_contours         = -0.1:0.02:0.1;
error_contours       = 0:0.01:0.1;
plot_sta             = 0;

load('usastates.mat');
statesdir = 'C:\Research\tstar\TDCode2D\tectonicprovinces\';

bnames={'Bound1'
        'Bound2'
        'Bound3'
        'Bound4'
        'Bound5'
        'Bound6'
        'Bound7'
        'Bound8'
        'Bound9'
        'Bound10'
        'Bound11'
        'Bound12'
        'Bound13'
        'Bound14'
        'Bound15'
        'Bound16'
        'Bound17'
        'Bound18'};
%%%%%%%
%set up

f = dir(loadname);

load( [ f(end).folder '\' f(end).name ])
a = whos('model_stats');

if isempty(a)
    
    allmodels = [];
    
    for k = 1:length(f)

        load([ f(k).folder '\' f(k).name ])
        m = models(:);

        allmodels = [ allmodels; m ];

    end
    
    centerLat = min(dataStruct.allLats)+(max(dataStruct.allLats)-min(dataStruct.allLats))/2;
    centerLon = min(dataStruct.allLons)+(max(dataStruct.allLons)-min(dataStruct.allLons))/2;

    origin               = [ centerLat centerLon ];           % array center [ lat lon ]
    mstruct              = defaultm('mercator');
    mstruct.origin       = [ origin TD_parameters.rotation ];%second number is rotation
    mstruct              = defaultm( mstruct );
    mstruct.scalefactor  = 6371;

    [xMat, yMat]    = meshgrid(dataStruct.xVec, dataStruct.yVec);
    [ model_stats ] = make_map_models( models, xMat, yMat, TD_parameters.sig_flag, -0.5:0.05:0.5, 0.66, TD_parameters.interp_style);

else
    
    load(loadname)
    
end
    
%%%%%%%
%tS map
figure(1)
%subplot(121)

m = model_stats.mean;
m(m < min(map_contours)) = min(map_contours);
m(m > max(map_contours)) = max(map_contours);

contourf(xMat, yMat, m, map_contours);
h = colorbar; h.Label.String = label;%h.Label.String = '\Deltat*';
hold on
plot(xMat(model_stats.std > masking_threshold), yMat(model_stats.std > masking_threshold), 'ws', 'MarkerFaceColor', 'w'); 

if plot_sta
   
    plot(dataStruct.dataX, dataStruct.dataY, 'ko')
    
end

xl = xlim;
yl = ylim;

if states_or_boundaries == 1

    for k = 1:length(usastates)

        [sx, sy] = mfwdtran(mstruct,usastates(k).Lat,usastates(k).Lon);
        plot(sx,sy,'k-','LineWidth',2)

    end

elseif states_or_boundaries == 2
    
    for knt=1:length(bnames)  

        B      = load([statesdir bnames{knt}]);
        [X, Y] = mfwdtran(mstruct,B(:,1),B(:,2));
        plot(X,Y,'k-','LineWidth',2);

    end

end
    
xlim(xl);
ylim(yl);

set(gca, 'XTick', []);
set(gca, 'YTick', []);

caxis([min(map_contours) max(map_contours)])

colormap(cmap)

%%%%%%
%std
figure(2)
%subplot(122)
m = model_stats.std;
m(m < min(error_contours)) = min(error_contours);
m(m > max(error_contours)) = max(error_contours);

contourf(xMat, yMat, m, error_contours);
h = colorbar; h.Label.String = [ '\sigma(' label ')' ];
hold on

if plot_sta
   
    plot(dataStruct.dataX, dataStruct.dataY, 'ko')
    
end

xl = xlim;
yl = ylim;

if states_or_boundaries == 1

    for k = 1:length(usastates)

        [sx, sy] = mfwdtran(mstruct,usastates(k).Lat,usastates(k).Lon);
        plot(sx,sy,'k-','LineWidth',2)

    end

elseif states_or_boundaries == 2
    
    for knt=1:length(bnames)  

        B      = load([statesdir bnames{knt}]);
        [X, Y] = mfwdtran(mstruct,B(:,1),B(:,2));
        plot(X,Y,'k-','LineWidth',2);

    end

end
    
xlim(xl);
ylim(yl);

set(gca, 'XTick', []);
set(gca, 'YTick', []);

colormap(cmap)

save(savename)
