clear
clc
close all

addpath('.\wfTools');

name = 'TAEmbayment_P_047';

data_to_load = [ '../Data/' name '.mat' ];
name_for_run = [ name '_results' ];%do not include .mat

tStarV = -0.25:0.01:1;%note that the final measurement is demeaned
pick_data     = 1;%if you are rerunning picked data set these to zero
repick_window = 1;
confidence_interval = 0.66;
xlimits = [ 605 620 ];%in seconds
data_cut_limits = [5 5];%in seconds
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

if pick_data

    load(data_to_load)

    chan = {Traces.channel};
    Traces(~strcmp(chan, 'BHZ')) = [];
    
    %need to remove instrument response first
    Traces = wfRemInstResp(Traces);
        
    lon = [Traces.longitude];
    [lon, ind] = sort(lon);
    Traces = Traces(ind);
    
    if ~isnan(filter_bounds)
        
        Traces = wfCosFilter( Traces, filter_bounds);

    end
            
    delete_vector = zeros(1, length(Traces));

    min_t = 1e9;
    
    %pad traces with 10 seconds of zeros in either direction
    for k = 1:length(Traces)

        Traces(k).data = [zeros(10*Traces(1).sampleRate,1); Traces(k).data; zeros(10*Traces(1).sampleRate,1)];
            
        min_t = min([ min_t length(Traces(k).data) ]);
        
    end
    
    for k = 1:length(Traces)

        Traces(k).data = Traces(k).data(1:min_t);
        
    end
    
    t = (1:length(Traces(1).data))/Traces(1).sampleRate;
    t_window = (t > xlimits(1)) & (t < xlimits(2));
    
    pre=data_cut_limits(1)*Traces(1).sampleRate;
    post=data_cut_limits(2)*Traces(1).sampleRate;

    figure(500)
    subplot(211)
    hold on
    %plot all of the traces
    for k = 1:length(Traces)

        if flip_waveforms
            
            Traces(k).data = Traces(k).data*-1;
            
        end
        
        %convert to displacement
        %Traces(k).data = cumsum(Traces(k).data);
        Traces(k).data = Traces(k).data/max(Traces(k).data(t_window));
                
        line_handles(k) = plot(t, Traces(k).data(1:length(t)) + k, 'k');
        
    end
    
    xlabel('Time, s')
    xlim(xlimits)
    subplot(212)
    
    t_ind = find(t>=xlimits(1) & t<=xlimits(2));
    
    for k = 1:length(Traces)

        subplot(212)
        plot(t, Traces(k).data, 'k')
        title(Traces(k).station)
        xlim(xlimits);
        ylim([-1 1])
        subplot(211)
        line_handles(k).LineWidth = 2;
        line_handles(k).Color     = 'r';
        shg

        subplot(212)       
        disp('Click below to delete, above to select for source, and on it to keep it but not use in the source');
        [ix, iy] = ginput(1);
        %just pick the max
        %[~, ix] = max(Traces(k).data(t_ind));
        
        ix = ix*Traces(k).sampleRate;
        
        if iy < -1

            delete_vector(k) = 1;
            line_handles(k).Visible = 'off';

            Traces_tmp(:,k) = Traces(k).data(1:pre+post+1);%gets deleted soon, just a dumy

        else
                
            Traces(k).data          = Traces(k).data(round(ix) - pre:round(ix) + post);
            Traces(k).data          = Traces(k).data.*tukeywin(length(Traces(k).data), 0.99);
            
            if iy > 1
                
                source_count                  = source_count + 1;
                
                if source_count == 1
                
                    source_traces              = Traces(k);
                    
                else
                   
                    source_traces(source_count) = Traces(k);
                    
                end
                
            end
            
            pick_time(k) = t(round(ix));

            line_handles(k).LineWidth = 0.5;
            line_handles(k).Color     = 'k';
                        
        end
        
        if any(isnan(Traces(k).data))
            
            delete_vector(k) = 1;
            
        end

    end
    
    Traces(logical(delete_vector))            = [];
    %get an estimate of the source
    %source = estimateSource_tmp(Traces, pre+1);
    source = estimateSource_tmp(source_traces);

    figure(100)
    
    plot(source);
    xlabel('Samples');
    title('Pick the window to fit over')
    [fitting_window, ~] = ginput(2);
    
    fitting_window = round(fitting_window);
    
    save([ name_for_run '.mat' ], 'Traces', 'source', 'fitting_window', 'pre', 'filter_bounds');

else

    load([ name_for_run '.mat' ]);

end

fitting_window = round(fitting_window);

if repick_window && ~pick_data
    
    figure(100)
    
    plot(source);
    xlabel('Samples');
    title('Pick the window to fit over')
    [fitting_window, ~] = ginput(2);
    
    fitting_window = round(fitting_window);

end
    
load coast
figure(1000)
hold on
%plot(states_lon, states_lat, 'k', 'LineWidth',1);
%plot(long, lat, 'k', 'LineWidth', 2);
xlim(lonlim); ylim(latlim)

figure(2000)
hold on
%plot(states_lon, states_lat, 'k', 'LineWidth', 1);
%plot(long, lat, 'k', 'LineWidth', 2);
xlim(lonlim); ylim(latlim)

fitting_window = sort(fitting_window);%you need this for indexing

for k = 1:length(Traces)
   
    disp([ 'On Station ' Traces(k).station]);
    
    source_struct = Traces(k);
    source_struct.data = source';
    source_struct.sampleCount = length(source);
        
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


