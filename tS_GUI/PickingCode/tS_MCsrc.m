clear
clc
close all

addpath('.\wfTools');

name = 'TAEmbayment_P_061';

data_to_load = [ '../Data/' name '.mat' ];
name_for_run = [ name '_results' ];%do not include .mat

tStarV = -0.25:0.01:1;%note that the final measurement is demeaned
pick_data     = 1;%if you are rerunning picked data set these to zero
repick_window = 1;
confidence_interval = 0.66;
xlimits = [ 580 620 ];%in seconds
data_cut_limits = [10 10];%in seconds
misfit_cut = 0.01;%set to nan is larger
latlim = [-90 90];
lonlim = [-180 180];

%flip waveforms? For mostly negative traces
flip_waveforms = 1;

insta_error_window = 20;%in samples
if_plot_bounds     = [ 0 2 ];
filter_bounds      = [ 0.01 0.02 1.5 3 ];

colorlimits = [-0.2  0.2];

source_count  = 0;

load('TAEmbayment_037_MC.mat');

[~, ind] = min([bestmodels.phi]);

source_struct = Traces(1);

wf = interp1(bestmodels(ind).tCell, bestmodels(ind).AmpCell, 1:length(Traces(1).data), 'spline');
fitting_window = [1 length(wf)];

for k = 1:length(Traces)
   
    disp([ 'On Station ' Traces(k).station]);
    
    source_struct = Traces(k);
    source_struct.data = wf';
    source_struct.sampleCount = length(wf);
        
    if ~isnan(filter_bounds)
        
        source_struct = wfCosFilter( source_struct, filter_bounds);

    end
    
    [ts_run(k), misfits] = tsGridSearch_jsb(Traces(k),source_struct,tStarV, pre, fitting_window, filter_bounds);
    t = (1:length(ts_run(k).data))/Traces(1).sampleRate;
            
end

figure(8)
hold on
title('Waveforms (black) against synthetics (red)');
xlabel('Time, s');
ylabel('Longitude of Station');

for k = 1:length(ts_run)

    t = (1:length(ts_run(k).data))/Traces(1).sampleRate - data_cut_limits(1);
        
    ts_vector(k) = ts_run(k).tStar;
    
    plot(t, ts_run(k).longitude + Traces(k).data(1:length(t))/max(Traces(k).data(1:length(t)))/15, 'k');
    plot(t, ts_run(k).longitude + ts_run(k).data/max(ts_run(k).data)/15, 'r');
    drawnow

end
 
ts_mean = nanmean(ts_vector);
ts_vector = ts_vector - ts_mean;

figure(1000)
hold on
figure(2000)


%demean all of the ts results?
for k = 1:length(Traces)
           
    figure(1000)
    
    if ~isnan(ts_vector(k))
    
        scatter(ts_run(k).longitude, ts_run(k).latitude, 30, ts_vector(k), 'filled');
        %scatter(ts_run(k).longitude, ts_run(k).latitude, 30, insta_freq_pick(k), 'filled');
        
    else
       
        scatter(ts_run(k).longitude, ts_run(k).latitude, 30, 'w', 'filled', 'MarkerEdgeColor', 'k');
        
    end
    
    figure(2000)
    scatter(ts_run(k).longitude, ts_run(k).latitude, 30, ts_run(k).misfit, 'filled');
    
    misfit(k) = ts_run(k).misfit;
    
end

figure(2000)
h = colorbar;
h.Label.String = 'Misfit';
title('Misfit');
caxis([0 misfit_cut])

figure(3000)
plot([ts_run.longitude], ts_vector, 'k.', 'LineStyle', 'none');
xlabel('Longitude of station');
ylabel('\Deltat*, s');
%line([-115.8 -115.8], [ -1 1 ], 'LineStyle', '--', 'LineWidth', 1, 'Color', 'k');
title('\Delta t* by longitude')
ylim([tStarV(1) tStarV(end)])

% figure(5000)
% plot(ts_vector, [ts_run.misfit], 'k.', 'LineStyle', 'none');
% xlabel('\Deltat*, s');
% ylabel('Misfit');
% title('\Delta t* against misfits')


save([name_for_run '_result.mat'])


