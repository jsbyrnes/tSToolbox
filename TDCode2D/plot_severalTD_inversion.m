%fit a smooth surface to all t* observations through an inversion scheme.
clear
close all
clc

load('AllRio_linear')

prefix = './AllRio_linearinversion*';

f = dir(prefix);

bm = [];

for k = 1:length(f)
   
    load(f(k).name)
    m = models(:);
    
    bm = [ bm; m ];
    
end

[xMat, yMat]    = meshgrid(dataStruct.xVec, dataStruct.yVec);
[ model_stats ] = make_map_models( bm, xMat, yMat, TD_parameters.sig_flag, -0.5:0.05:0.5, 0.66, TD_parameters.interp_style);

m = model_stats.mean;
m(m<-0.1) = -0.1;

figure(1)
subplot(211)
contourf(xMat, yMat, m, -0.1:0.02:0.1);
h = colorbar; h.Label.String = '\Deltat*';
hold on
%plot(dataStruct.dataX, dataStruct.dataY, 'rs');
mask_x = xMat(model_stats.std > 0.15);
mask_y = yMat(model_stats.std > 0.15);

plot(mask_x, mask_y, 'o', 'MarkerFaceColor', [ 1 1 1 ], 'MarkerEdgeColor', [ 1 1 1 ])

subplot(212)
contourf(xMat, yMat, model_stats.std, 0.00:0.02:0.2);
h = colorbar; h.Label.String = '\sigma(\Deltat*)';
hold on
plot(dataStruct.dataX, dataStruct.dataY, 'rs');

[X, P] = meshgrid(xMat(1, :), -0.5:0.05:0.5);

figure
contourf(X', P', squeeze(model_stats.pdf_set(3,:,:)), 10)
colorbar
colormap(flipud(pink))
