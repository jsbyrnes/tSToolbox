function varargout = sourceTracesSelection(varargin)
% SOURCETRACESSELECTION MATLAB code for sourceTracesSelection.fig
%      SOURCETRACESSELECTION, by itself, creates a new SOURCETRACESSELECTION or raises the existing
%      singleton*.
%
%      H = SOURCETRACESSELECTION returns the handle to a new SOURCETRACESSELECTION or the handle to
%      the existing singleton*.
%
%      SOURCETRACESSELECTION('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SOURCETRACESSELECTION.M with the given input arguments.
%
%      SOURCETRACESSELECTION('Property','Value',...) creates a new SOURCETRACESSELECTION or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before sourceTracesSelection_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to sourceTracesSelection_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help sourceTracesSelection

% Last Modified by GUIDE v2.5 06-Aug-2019 12:04:56

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @sourceTracesSelection_OpeningFcn, ...
                   'gui_OutputFcn',  @sourceTracesSelection_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end

addpath('./wfTools/')
addpath('./PickingCode/')

% End initialization code - DO NOT EDIT


% --- Executes just before sourceTracesSelection is made visible.
function sourceTracesSelection_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to sourceTracesSelection (see VARARGIN)

% Choose default command line output for sourceTracesSelection
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes sourceTracesSelection wait for user response (see UIRESUME)
% uiwait(handles.figure1);

dataLoaded = 0;
setappdata(gcf, 'dataLoaded', dataLoaded);
setappdata(gcf, 'plotmode', 1);
addpath(genpath('./wfTools/'))
addpath(genpath('./PickingCode/'))

addpath('../irisFetch/')
javaaddpath('../irisFetch/IRIS-WS-2.0.18.jar')

% set the keypress_fcn for the figure
set(gcf,'WindowKeyPressFcn',@chngTr)

% --- Outputs from this function are returned to the command line.
function varargout = sourceTracesSelection_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

% --- Executes on button press in done_btn.
function done_btn_Callback(hObject, eventdata, handles)
% hObject    handle to done_btn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    
    %keep only the traces that were not culled out
    Traces=getappdata(gcf,'Traces');
    %useVec=getappdata(gcf,'useVec');
    %Traces=Traces(useVec);
    fw_start       = getappdata(gcf,'fw_start');
    fitting_window = getappdata(gcf,'fitting_window');
    filter_bounds  = getappdata(gcf,'filter_bounds');
    t              = getappdata(gcf,'t');
    xrange         = getappdata(gcf,'xrange');
    midpoint       = getappdata(gcf,'midpoint');

    %plot source and pick fitting window
    forSrc=getappdata(gcf,'forSrc');
    srcTr=forSrc{2};
    MsrcTr=mean(srcTr,2); %stack of selected traces

    source=MsrcTr; %just renaming

    %cut the data and the source
    for i = 1:length(Traces)

        Traces(i).data        = Traces(i).data/rms(Traces(i).data(t > fw_start & t < fitting_window(2)));
        Traces(i).data        = Traces(i).data(t > (midpoint - xrange) & t < (midpoint + xrange));
        Traces(i).sampleCount = length(Traces(i).data);

    end

    source         = source((t > (midpoint - xrange) & t < (midpoint + xrange)));
    [~, ind]       = min(abs(t - (midpoint - xrange)));
    fw_start       = round((fw_start - t(ind))*Traces(1).sampleRate);
    fitting_window = round((fitting_window - t(ind))*Traces(1).sampleRate);

    if fw_start(1) < 1
        
        fw_start(1) = 1;
        
    end
    
    if fitting_window(1) < 1
        
        fitting_window(1) = 1;
        
    end
    
    if fitting_window(2) > length(t)
       
        fitting_window(2) = length(t);
        
    end
    
    tStarV           = -0.1:0.05:2;%note that the final measurement is demeaned
    
    fw_ending = (fitting_window(1) + 1):fitting_window(2);
    
    for k = 1:length(fw_ending)

        fprintf([ 'On window ' num2str(k) ' of ' num2str(length(fw_ending)) ' candidate windows\n' ]);
        
        for kk = 1:length(Traces)

            source_struct = Traces(kk);
            source_struct.data = source;
            source_struct.sampleCount = length(source);

            %if ~isnan(filter_bounds)

            %   source_struct = wfButterworth( source_struct, filter_bounds);

            %end

            [ts_run(k, kk), ~] = tsGridSearch_interp(Traces(kk), source_struct, tStarV, [ fw_start fw_ending(k) ], filter_bounds);
            %[ts_run(k, kk), ~] = tsGridSearch_jsb(Traces(kk), source_struct, tStarV, [ fw_start fw_ending(k) ], filter_bounds);

        end
    
    end
    
    setappdata(gcf, 'ts_run', ts_run);
    setappdata(gcf, 'Traces', Traces);
    setappdata(gcf, 'fw_end', fitting_window(2));
    handles.FWslide.Max         = diff(fitting_window);
    handles.FWslide.SliderStep  = [ 1/handles.FWslide.Max 3/handles.FWslide.Max ];
    handles.FWslide.Value       = handles.FWslide.Max;
    
    lh = getappdata(gcf, 'lineHandles');
    delete(lh)

    %resize the source traces for plotting beforfe resizing t
    srcTr = srcTr((t > (midpoint - xrange) & t < (midpoint + xrange)), :);

    t = (1:length(Traces(1).data))/Traces(1).sampleRate;

    axes(handles.allWf_ax)
    hold on
    %plot all of the traces
    for k = 1:length(Traces)

        Traces(k).data = Traces(k).data/max(abs(Traces(k).data(fw_start:fitting_window(2))));

        line_handles(k) = plot(t, Traces(k).data + k, 'k-', 'lineWidth', 1);

    end

    xlabel('Time, s')
    xlim([1 max(t)])
    ylim([-1 (length(line_handles)+1)])

    setappdata(gcf,'lineHandles',line_handles);
    setappdata(gcf,'t',t);

    fwh            = getappdata(gcf, 'fwh');

    delete(fwh);

    forSrc{2} = srcTr;
    setappdata(gcf,'forSrc',forSrc);
    fw_start       = fw_start/Traces(1).sampleRate;
    fitting_window = fitting_window/Traces(1).sampleRate;
    setappdata(gcf, 'fw_start', fw_start);
    setappdata(gcf, 'fitting_window', fitting_window);

    %update plotting variables
    midpoint = mean(t);
    xrange   = max(t)/2;

    setappdata(gcf, 'midpoint', midpoint);
    setappdata(gcf, 'xrange', xrange);

    axes(handles.currWf_ax)
    ylimits = ylim;
    stfwh   = plot( [ fw_start fw_start ], [ ylimits(1) ylimits(2) ], 'k--', 'LineWidth', 0.5);
    fwh(1)  = plot( [ fitting_window(1) fitting_window(1) ], [ ylimits(1) ylimits(2) ], 'k--', 'LineWidth', 0.5);
    fwh(2)  = plot( [ fitting_window(2) fitting_window(2) ], [ ylimits(1) ylimits(2) ], 'k--', 'LineWidth', 0.5);
    setappdata(gcf, 'stfwh', stfwh);
    setappdata(gcf, 'fwh', fwh);

    handles.QCButton.Enable            = 'on';
    handles.QCButton.Value             = true;
    handles.DataSelectionButton.Value  = false;
        
    QCButton_Callback(handles.QCButton, [], handles)
    DataSelectionButton_Callback(handles.DataSelectionButton, [], handles)
    
    pltSrcTrs
    plotMap
    
function chngTr(figH,key)

dataLoaded = getappdata(gcf, 'dataLoaded');

if dataLoaded
    
    handles  = guihandles;
    inSrcVec = getappdata(gcf,'inSrcVec');
    plotmode = getappdata(gcf,'plotmode');

    %the following is so that traces included in source estimation are
    %highlighted
    lw=ones(size(inSrcVec));

    if plotmode == 1
        
        lw(inSrcVec)=4;
        lc = 'b';
        
    else
       
        lc = 'k';
        
    end
    if isempty(key.Modifier)

        if strcmp(key.Key,'uparrow')
            cwf=getappdata(gcf,'currentWf');
            lh=getappdata(gcf,'lineHandles');
            nTraces=length(lh);

            if cwf<nTraces

                set(lh(cwf),'color','k','lineWidth',lw(cwf))
                cwf=cwf+1;
                set(lh(cwf),'color',lc,'lineWidth',4)
                setappdata(gcf,'currentWf',cwf)
                
                pltSrcTrs
                
            end

        elseif strcmp(key.Key,'downarrow')
            
            cwf=getappdata(gcf,'currentWf');
            lh=getappdata(gcf,'lineHandles');

            if cwf>1

                set(lh(cwf),'color','k','lineWidth',lw(cwf))
                cwf=cwf-1;
                set(lh(cwf),'color',lc,'lineWidth',4)
                setappdata(gcf,'currentWf',cwf)

                pltSrcTrs

            end

        elseif strcmp(key.Key,'leftarrow')

            cwf=getappdata(gcf,'currentWf');
            lh=getappdata(gcf,'lineHandles');

            Traces = getappdata(gcf,'Traces');
            Traces(cwf).data = circshift(Traces(cwf).data, -1);        
            lh(cwf).YData = Traces(cwf).data + cwf;
            setappdata(gcf, 'Traces', Traces);
                       
            add_to_source%remove the old
            add_to_source%add the shifted
                        
            pltSrcTrs

        elseif strcmp(key.Key,'rightarrow')

            cwf=getappdata(gcf,'currentWf');
            lh=getappdata(gcf,'lineHandles');

            Traces = getappdata(gcf,'Traces');
            Traces(cwf).data = circshift(Traces(cwf).data, 1);        
            lh(cwf).YData = Traces(cwf).data + cwf;
            setappdata(gcf, 'Traces', Traces);

            add_to_source%remove the old
            add_to_source%add the shifted
                        
            pltSrcTrs
            
        elseif strcmp(key.Key,'backspace')

            useFunc%flips if in or out of useVec
            pltSrcTrs
            
        elseif strcmp(key.Key,'space')

            add_to_source%flips if in or out of source
            pltSrcTrs
            
        end

    elseif ~isempty(key.Modifier)

        if strcmp(key.Modifier{1}, 'shift')%adjust the position of the window

            axes(handles.currWf_ax)
            midpoint = getappdata(gcf, 'midpoint');
            xrange   = getappdata(gcf, 'xrange');
            t        = getappdata(gcf, 't');

            if strcmp(key.Key,'leftarrow')

                midpoint = midpoint - 1;

            elseif strcmp(key.Key,'rightarrow')

                midpoint = midpoint + 1;

            elseif strcmp(key.Key,'uparrow')

                xrange = xrange + xrange*0.05;

            elseif strcmp(key.Key,'downarrow')

                xrange = xrange - xrange*0.05;

            end

            minval = max([ 0 (midpoint - xrange) ]);
            maxval = min([ (midpoint + xrange) max(t) ]);
            xlim([ minval maxval ]);

            axes(handles.allWf_ax)
            xlim([ minval maxval ]);

            setappdata(gcf,'midpoint', midpoint); %these are the traces that will go into the source estimate
            setappdata(gcf,'xrange',   xrange); %these are the traces that will go into the source estimate
            
        end
        
    end
    
    plotMap
    
end

function pltSrcTrs

handles  = guihandles;
plotmode = getappdata(gcf, 'plotmode');
useVec   = logical(getappdata(gcf, 'useVec'));
Traces   = getappdata(gcf, 'Traces');
cwf      = getappdata(gcf, 'currentWf');
t        = getappdata(gcf, 't');
forSrc   = getappdata(gcf, 'forSrc');
stfw     = getappdata(gcf, 'fw_start');
fw       = getappdata(gcf, 'fitting_window');


if useVec(cwf) == 1
    
    ls = '-';
    
else
    
    ls = '--';
    
end

axes(handles.currWf_ax)

if plotmode == 1
    
    srcTr  = forSrc{2};
    
    %delete stuff if something already there
    h=findobj('Tag','srcTr');
    if ~isempty(h); delete(h); end
    h=findobj('Tag','srcTrMany');
    if ~isempty(h); delete(h); end
    h=findobj('Tag','synTr');
    if ~isempty(h); delete(h); end
    h=findobj('Tag','theTrace');
    if ~isempty(h); delete(h); end
    
    stfwhtmp    = getappdata(gcf, 'stfwh');
    stfwhalltmp = getappdata(gcf, 'fwh');
    
    if ~isempty(stfwhtmp)    delete(stfwhtmp);    end
    if ~isempty(stfwhalltmp) delete(stfwhalltmp); end
    
    if ~isempty(srcTr)
                        
        for k = 1:length(forSrc{1})
        
            srcTr(:, k) = srcTr(:, k)...
                /rms(srcTr(t > stfw & t < fw(2), k));
            
        end
        
        h=plot(t, srcTr,'-','color',[.6 .6 .6]);
        set(h,'Tag','srcTrMany');
        
        MsrcTr=mean(srcTr,2);
        h=plot(t, MsrcTr,'-','color','k','linewidth',1.5);
        set(h,'Tag','srcTr');
        
        delete(findobj('Tag','theTrace')); %this just gets rid of the previous one
        h=plot(t, Traces(cwf).data/rms(Traces(cwf).data(t > stfw & t < fw(2))),'b', 'LineStyle', ls, 'linewidth',1.5);
        set(h,'Tag','theTrace')
 
%         h = legend('All source traces', 'Current source estimate', 'Current data', 'fwStart', 'fwEnd', 'fwMax');
        h = legend('All source traces', 'Current source estimate', 'Current data');
        set(h, 'Tag', 'srcLegend');
        set(h, 'Location', 'northwest');
        
        ylimits = ylim;
        stfwh = plot( [ stfw stfw ], [ ylimits(1) ylimits(2) ], 'k--', 'LineWidth', 1.5);
        setappdata(gcf, 'stfwh', stfwh);

        fwh(1) = plot( [ fw(1) fw(1) ], [ ylimits(1) ylimits(2) ], 'k--', 'LineWidth', 0.5);
        fwh(2) = plot( [ fw(2) fw(2) ], [ ylimits(1) ylimits(2) ], 'k--', 'LineWidth', 0.5);
        setappdata(gcf, 'fwh', fwh);
        
    else
                
        delete(findobj('Tag','theTrace')); %this just gets rid of the previous one
        h=plot(t, Traces(cwf).data/rms(Traces(cwf).data(t > stfw & t < fw(2))),'b', 'LineStyle', ls, 'linewidth',1.5);
        set(h,'Tag','theTrace')
        
        ylimits = ylim;
        stfwh = plot( [ stfw stfw ], [ ylimits(1) ylimits(2) ], 'k--', 'LineWidth', 1.5);
        setappdata(gcf, 'stfwh', stfwh);

        fwh(1) = plot( [ fw(1) fw(1) ], [ ylimits(1) ylimits(2) ], 'k--', 'LineWidth', 0.5);
        fwh(2) = plot( [ fw(2) fw(2) ], [ ylimits(1) ylimits(2) ], 'k--', 'LineWidth', 0.5);
        setappdata(gcf, 'fwh', fwh);
        
    end
    
elseif plotmode == 2
    
    ts_run = getappdata(gcf, 'ts_run');
    fw_ind = round(handles.FWslide.Value);
    
    [nfw, ~] = size(ts_run);
    
    %delete stuff if something already there
    h=findobj('Tag','synTr');
    if ~isempty(h); delete(h); end
    h=findobj('Tag','srcTr');
    if ~isempty(h); delete(h); end
    h=findobj('Tag','srcTrMany');
    if ~isempty(h); delete(h); end
    h=findobj('Tag','theTrace');
    if ~isempty(h); delete(h); end
    h=findobj('Tag','srcLegend');
    if ~isempty(h); delete(h); end
    
    stfwh            = getappdata(gcf, 'stfwh');
    fwh            = getappdata(gcf, 'fwh');
    if ~isempty(stfwh); delete(stfwh); end
    if ~isempty(fwh); delete(fwh); end
    
    if ~isempty(ts_run)
        
        h=plot(t, ts_run(fw_ind, cwf).data/rms(ts_run(fw_ind, cwf).data(t > stfw & t < fw(2))),'r', 'LineStyle', ls, 'linewidth', 1.5);
        set(h,'Tag','synTr');
        
        h=plot(t, Traces(cwf).data/rms(Traces(cwf).data(t > stfw & t < fw(2))),'-','color','k', 'LineStyle', ls, 'linewidth',1.5);
        set(h,'Tag','srcTr');
        
    end
    
%     h = legend('Data', [ 'Synthetic: \Deltat* = ' num2str(ts_run(fw_ind, cwf).tStar_WF, 3)], 'fwStart', 'fwEnd', 'fwMax');
    h = legend('Data', [ 'Synthetic: \Deltat* = ' num2str(ts_run(fw_ind, cwf).tStar_WF, 3)]);
    set(h, 'Tag', 'srcLegend');
    set(h, 'Location', 'northwest');
        
    fw = getappdata(gcf, 'fitting_window');
    %fwh(2).XData   = [ (fw(1) + (handles.FWslide.Value/Traces(1).sampleRate)) (fw(1) + (handles.FWslide.Value/Traces(1).sampleRate)) ];

    ylimits = ylim;
    stfwh = plot( [ stfw stfw ], [ ylimits(1) ylimits(2) ], 'k--', 'LineWidth', 1.5);
    setappdata(gcf, 'stfwh', stfwh);
    fwh(1) = plot( [ fw(1) fw(1) ], [ ylimits(1) ylimits(2) ], 'k--', 'LineWidth', 0.5);
    fwh(2) = plot( [ (fw(1) + (handles.FWslide.Value/Traces(1).sampleRate)) (fw(1) + (handles.FWslide.Value/Traces(1).sampleRate)) ], [ ylimits(1) ylimits(2) ], 'k--', 'LineWidth', 0.5);
    fwh(1).Visible = 'off';
    setappdata(gcf, 'fwh', fwh);
    
    axes(handles.tSFWplot)
    hold on
    htSFW  = getappdata(gcf, 'htSFW');
    htSFWl = getappdata(gcf, 'htSFWl');
    if ~isempty(htSFW);  delete(htSFW);  end
    if ~isempty(htSFWl); delete(htSFWl);  end

    xax = linspace(fw(1), fw(2), nfw);
    htSFW = plot(xax, [ts_run(:, cwf).tStar_WF], 'ko');
    htSFWl = plot([xax(fw_ind) xax(fw_ind) ], ylim, 'k--');
    xlabel('Fitting window ending, s');
    ylabel('\Deltat*, s');
    setappdata(gcf, 'htSFW',  htSFW);
    setappdata(gcf, 'htSFWl', htSFWl);
    
end

axes(handles.currWf_ax)

midpoint = getappdata(gcf, 'midpoint');
xrange   = getappdata(gcf, 'xrange');
xlim([ (midpoint - xrange) (midpoint + xrange) ]);
%ylim([-1 1]);

delete(findall(gcf,'type','annotation'));
textLoc=[0.55 0.5 0.3 0.3];
annotation('textbox', textLoc, 'String', [Traces(cwf).network,'.',Traces(cwf).station],'FitBoxToText','on'); 

set(gca, 'YTick', []);
   
function plotMap

    handles  = guihandles;

    useVec          = logical(getappdata(gcf, 'useVec'));
    ts_run          = getappdata(gcf, 'ts_run');
    cwf             = getappdata(gcf, 'currentWf');
        
    fw_ind = round(handles.FWslide.Value);

    if isempty(ts_run)
        
        ts_run = getappdata(gcf, 'Traces');
        fw_ind = 1;
        
        for k = 1:length(ts_run)
           
            ts_run(k).tStar_WF = 0;
            
        end
        
    end
    
    axes(handles.map_ax)
    hold on

    hMapQC  = getappdata(gcf, 'hMapQC');
    hMapNQC = getappdata(gcf, 'hMapNQC');
    hMark   = getappdata(gcf, 'hMark');
    hC      = getappdata(gcf, 'hC');
    if ~isempty(hMapQC);    delete(hMapQC);  end
    if ~isempty(hMapNQC);   delete(hMapNQC); end
    if ~isempty(hMark);     delete(hMark);   end
    if ~isempty(hC);        delete(hC);   end
    
    load('CMfine')
    hMapQC  = scatter( [ts_run(fw_ind, useVec).longitude],  [ts_run(fw_ind, useVec).latitude],  20, [ts_run(fw_ind, useVec).tStar_WF],  'filled');
    hMapNQC = scatter( [ts_run(fw_ind, ~useVec).longitude], [ts_run(fw_ind, ~useVec).latitude], 20, [ts_run(fw_ind, ~useVec).tStar_WF], 'filled');
    %colormap(cm); hC = colorbar; hC.Location = 'westoutside'; hC.Label.String = '\Deltat*, s';
    %caxis([ min([ts_run(fw_ind, useVec).tStar_WF] - 0.1) max([ts_run(fw_ind, useVec).tStar_WF] + 0.1)]);
    
    colormap(cm); 

    % edited by Zhao, fix the problem of multi colorbars when press
    % backspace (when changing the file path)
    C = findall(gcf,'type','ColorBar');
    if isempty(C)
        hC = colorbar; hC.Location = 'westoutside'; hC.Label.String = '\Deltat*, s';
        if(any([ts_run(fw_ind, useVec).tStar_WF])) % check if the field exists
            caxis([ min([ts_run(fw_ind, useVec).tStar_WF] - 0.1) max([ts_run(fw_ind, useVec).tStar_WF] + 0.1)]);
        else
            caxis([-0.1 0.1]);
        end
    end
    
    xlabel(['Longitude, ' char(176)]);
    ylabel(['Latitude, ' char(176)]);
    
    hMapQC.MarkerEdgeColor  = 'k';
    hMapNQC.MarkerFaceColor = 'k';
    
    hMapQC.ButtonDownFcn  = @ClickMap;
    hMapNQC.ButtonDownFcn = @ClickMap;
    
    hMark = plot(ts_run(1, cwf).longitude, ts_run(1, cwf).latitude, 'k*', 'MarkerSize', 10);

    setappdata(gcf,'hMapQC',  hMapQC);
    setappdata(gcf,'hMapNQC', hMapNQC);
    setappdata(gcf,'hMark',   hMark);
    setappdata(gcf,'hC',      hC);
    
function add_to_source
      
    inSrcVec = getappdata(gcf,'inSrcVec');
    cwf      = getappdata(gcf,'currentWf');
    forSrc   = getappdata(gcf,'forSrc');
    srcIx    = forSrc{1}; srcTr=forSrc{2};

    if ~inSrcVec(cwf)
        %add this wf number to the list
       srcIx=[srcIx; cwf]; 
       %add this wf to the set
       Traces=getappdata(gcf,'Traces');
       currWf=Traces(cwf).data;%was 75/125
       srcTr=[srcTr currWf];
       forSrc={srcIx srcTr};
       setappdata(gcf,'forSrc',forSrc)
       inSrcVec(cwf)=true;

    else

       %find the column number in srcTr where the current trace is stored and
       %get rid of it, then put the things back into the structure and into
       %appdata
       inSrcVec(cwf)=false; 

       ixx=find(srcIx==cwf);

       %remove the index and trace from the source lists
       srcIx(ixx)=[];
       srcTr(:,ixx)=[];

       forSrc={srcIx srcTr};
       setappdata(gcf,'forSrc',forSrc)

    end

    setappdata(gcf,'inSrcVec',inSrcVec);
    
function useFunc

    useVec      = logical(getappdata(gcf,'useVec'));
    cwf         = getappdata(gcf,'currentWf');
    lineHandles = getappdata(gcf,'lineHandles');

    if ~useVec(cwf)

       useVec(cwf)=true;
       set(lineHandles(cwf),'lineStyle','-')

    else

       useVec(cwf)=false; 
       set(lineHandles(cwf),'lineStyle','--')

    end

    setappdata(gcf,'useVec',useVec);

% --- Executes on button press in SaveDataButton.
function SaveDataButton_Callback(hObject, eventdata, handles)
% hObject    handle to SaveDataButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
ts_run_all   = getappdata(gcf, 'ts_run');
Traces       = getappdata(gcf, 'Traces');
useVec       = logical(getappdata(gcf, 'useVec'));
inSrcVec     = getappdata(gcf, 'inSrcVec');
fw_start     = getappdata(gcf, 'fw_start');
fw           = getappdata(gcf, 'fitting_window');

fw_ind = round(handles.FWslide.Value);
ts_run = squeeze(ts_run_all(fw_ind, :));

handles = guihandles;

for k = 1:length(Traces)
   
    Traces(k).QC    = useVec(k);
    Traces(k).inSrc = inSrcVec(k);
    
end

savename = handles.DataLoadBox.String;

s = strsplit(savename, '.');
if strcmp('mat', s(end))
      
    savename = savename(1:end-4);%remove the extension
        
end

save([savename 'Measurement.mat'], 'ts_run_all', 'Traces', 'ts_run', 'fw_start', 'fw', 'fw_ind');

function DataLoadBox_Callback(hObject, eventdata, handles)
% hObject    handle to DataLoadBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of DataLoadBox as text
%        str2double(get(hObject,'String')) returns contents of DataLoadBox as a double
datapath = get(hObject, 'String');
setappdata(gcf, 'datapath', datapath);
handles = guihandles;
handles.LoadDataButton.String = 'Load Data';
setappdata(gcf, 'dataLoaded', 0);

% --- Executes during object creation, after setting all properties.
function DataLoadBox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to DataLoadBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in LoadDataButton.
function LoadDataButton_Callback(hObject, eventdata, handles)
% hObject    handle to LoadDataButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% first some parameters and things, not sure I need all this.

handles  = guihandles;
name     = getappdata(gcf, 'datapath');
dl       = getappdata(gcf, 'dataLoaded');

if ~dl

    s = strsplit(name, '.');
    if ~strcmp('mat', s(end))

        name = [ name '.mat' ];%add the extension

    end
    
    xrange        = 30;    
    %get the data and plot it all

    if exist(name, 'file')

        load(name)
        
        handles.LoadDataButton.String = 'Loading...';
        drawnow

        chan = {Traces.channel};
        Traces(~((strcmp(chan, 'BHZ') | (strcmp(chan, 'HHZ'))))) = [];
        
        lon = [Traces.longitude];
        [~, sind] = sort(lon);
        Traces = Traces(sind);
        
        % remove instrument response only when it has not been removed
        % during fetching
        if(any([Traces.instrument]))
            Traces = wfRemInstResp(Traces); 
        end

        filter_bounds = [ str2num(handles.LowHz.String) str2num(handles.HighHz.String) ];
    
        filter_bounds(filter_bounds==0) = NaN;
        
        for k = 1:length(Traces)

            Traces(k).data = cumsum(Traces(k).data);
            
        end
        
        if ~isnan(filter_bounds)

           Traces = wfButterworth( Traces, filter_bounds);

        end

        for k = 1:length(Traces)

            Traces(k).data = diff(Traces(k).data);
            Traces(k).data = Traces(k).data - mean(Traces(k).data);
            Traces(k).data = Traces(k).data/rms(Traces(k).data);
            
        end
        
        % equalize trace length first
        minlen = 1e100;
        for k = 1:length(Traces)

            minlen = min([minlen length(Traces(k).data)]);

        end

        for k = 1:length(Traces)

            Traces(k).data = Traces(k).data(1:minlen);
            Traces(k).sampleCount = minlen;

        end

        xd             = minlen/Traces(1).sampleRate;
        midpoint       = (xd/2);
        xlimits        = [ (midpoint - xrange) ((midpoint) + xrange) ];
        fitting_window = [ (midpoint - xrange/2) (midpoint + xrange/2) ];
        fw_start       = (midpoint - xrange);
        
        t = (1:minlen)/Traces(1).sampleRate;
        
        lhtmp       = getappdata(gcf, 'lineHandles');
        stfwhtmp    = getappdata(gcf, 'stfwh');
        stfwhalltmp = getappdata(gcf, 'stfwh_all');
        fwhtmp      = getappdata(gcf, 'fwh');
        fwhalltmp   = getappdata(gcf, 'fwh_all');
        
        if ~isempty(lhtmp);       delete(lhtmp);       end
        if ~isempty(stfwhtmp);    delete(stfwhtmp);    end
        if ~isempty(stfwhalltmp); delete(stfwhalltmp); end
        if ~isempty(fwhtmp);      delete(fwhtmp);      end
        if ~isempty(fwhalltmp);   delete(fwhalltmp);   end
        
        % set "use" to yes for all of them unless they are out of range
        useVec = true(1, length(Traces));
        degVec = [ str2num(handles.MinDeg.String) str2num(handles.MaxDeg.String) ];
        
        for k = 1:length(Traces)
           
            [arc, ~] = distance(eventData.PreferredLatitude, eventData.PreferredLongitude, ...
                Traces(k).latitude, Traces(k).longitude);
            
            if arc < degVec(1) || arc > degVec(2)
                
                useVec(k) = 0;
                
            end
            
        end
        
        axes(handles.allWf_ax)
        hold on
        %plot all of the traces
        for k = 1:length(Traces)

            Traces(k).data = Traces(k).data/max(abs(Traces(k).data...
                (t > fw_start & t < fitting_window(2))));

            %variable name is changed because when you delete above, you
            %get an empty array. When you assign to the empty array, the
            %handle is a double, not an object. Seems a bug in matlab?
            
            if useVec(k)
            
                lh(k) = plot(t, Traces(k).data + k, 'k-', 'lineWidth', 1);
                                
            else
                
                lh(k) = plot(t, Traces(k).data + k, 'k--', 'lineWidth', 1);
                
            end

            lh(k).ButtonDownFcn = @ClickCWF;
            
        end

        xlabel('Time, s')
        xlim(xlimits)
        ylim([-1 (length(lh)+1)])
        
        % set "inSrc" to no for all of them
        inSrcVec = false(1, length(Traces));

        % plot the first trace in the currWf axis
        setappdata(gcf,'currentWf',     1);
        setappdata(gcf,'nTraces',       length(Traces));
        setappdata(gcf,'lineHandles',   lh);
        setappdata(gcf,'Traces',        Traces);
        setappdata(gcf,'useVec',        useVec);
        setappdata(gcf,'inSrcVec',      inSrcVec);
        setappdata(gcf,'midpoint',      midpoint); %these are the traces that will go into the source estimate
        setappdata(gcf,'xrange',        xrange); %these are the traces that will go into the source estimate
        setappdata(gcf,'t',             t);
        setappdata(gcf,'forSrc',        {[] []}); %these are the traces that will go into the source estimate
        setappdata(gcf,'ts_run',        []); %clear any preexisting_run
        setappdata(gcf,'level',         0); %clear any preexisting_run
        setappdata(gcf,'plotmode',      1);
        setappdata(gcf,'stfwh',         []);
        setappdata(gcf,'stfwh_all',     []);
        setappdata(gcf,'fwh',           []);
        setappdata(gcf,'fwh_all',       []);
        setappdata(gcf,'filter_bounds', filter_bounds);

        set(lh(1),'color','b','linewidth',4)

        axes(handles.currWf_ax); hold on
        h=plot(t, Traces(1).data/max(abs(Traces(k).data...
                (t > fw_start & t < fitting_window(2)))),'b-','linewidth',2);
        set(h,'Tag','theTrace')
        
        ylimits = ylim;
        stfwh  = plot( [ fw_start fw_start ], [ ylimits(1) ylimits(2) ], 'k--', 'LineWidth', 1.5);
        fwh(1) = plot( [ fitting_window(1) fitting_window(1) ], [ ylimits(1) ylimits(2) ], 'k--', 'LineWidth', 0.5);
        fwh(2) = plot( [ fitting_window(2) fitting_window(2) ], [ ylimits(1) ylimits(2) ], 'k--', 'LineWidth', 0.5);
        setappdata(gcf, 'stfwh', stfwh);
        setappdata(gcf, 'fwh', fwh);
        setappdata(gcf, 'fw_start', fw_start);
        setappdata(gcf, 'fitting_window', fitting_window);
        xlim(xlimits)
        setappdata(gcf, 'dataLoaded', 1);

        axes(handles.allWf_ax)
        ylimits = ylim;
        stfwh_all  = plot( [ fw_start fw_start ], [ ylimits(1) ylimits(2) ], 'k--', 'LineWidth', 1.5);
        fwh_all(1) = plot( [ fitting_window(1) fitting_window(1) ], [ ylimits(1) ylimits(2) ], 'k--', 'LineWidth', 0.5);
        fwh_all(2) = plot( [ fitting_window(2) fitting_window(2) ], [ ylimits(1) ylimits(2) ], 'k--', 'LineWidth', 0.5);
        setappdata(gcf, 'stfwh_all', stfwh_all);
        setappdata(gcf, 'fwh_all', fwh_all);

        handles.QCButton.Enable            = 'off';
        handles.QCButton.Value             = false;
        handles.DataSelectionButton.Value  = true;

        handles.LoadDataButton.String = 'Data Loaded';
        drawnow

        pltSrcTrs  
        plotMap
        
        handles.DiffButton.String = 'Acceleration';
        handles.DiffButton.Enable = 'on';
        handles.IntButton.String  = 'Displacement';
        handles.IntButton.Enable  = 'on';
        
    else

        handles.LoadDataButton.String = 'File not found';
        drawnow

    end
    
end
       
% --- Executes on button press in QCButton.
function QCButton_Callback(hObject, eventdata, handles)
% hObject    handle to QCButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if get(hObject, 'Value') == true
    
    setappdata(gcf, 'plotmode', 2);
    
    handles.DataSelectionButton.Value = false;

    cwf      = getappdata(gcf,'currentWf');
    lh       = getappdata(gcf,'lineHandles');

    useVec   = logical(getappdata(gcf,'useVec'));
        
    for k = 1:length(lh)

        if useVec(k)
        
            set(lh(k),'color','k', 'LineWidth', 1)
            
        else
           
            set(lh(k),'color','k', 'LineWidth', 1, 'LineStyle', '--')
            
        end

    end

    set(lh(cwf),'LineWidth',4)
    
    handles.RemoveButon.Enable = false;
   
elseif get(hObject, 'Value') == false
    
    handles.DataSelectionButton.Value = true;
    
end

pltSrcTrs

% Hint: get(hObject,'Value') returns toggle state of QCButton

% --- Executes on button press in DataSelectionButton.
function DataSelectionButton_Callback(hObject, eventdata, handles)
% hObject    handle to DataSelectionButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if get(hObject, 'Value') == true
    
    setappdata(gcf, 'plotmode', 1);
    handles.QCButton.Value = false;
    
    cwf      = getappdata(gcf,'currentWf');
    lh       = getappdata(gcf,'lineHandles');
    inSrcVec = getappdata(gcf,'inSrcVec');
    useVec   = logical(getappdata(gcf,'useVec'));

    lw = ones(size(lh));
    lw(inSrcVec) = 4;

    for k = 1:length(lh)

        if useVec(k)
        
            set(lh(k),'color','k','lineWidth',lw(k))
        
        else
            
            set(lh(k),'color','k','lineWidth',lw(k), 'LineStyle', '--')
            
        end

    end

    set(lh(cwf),'color','b')
    
    handles.RemoveButon.Enable = true;
    
elseif get(hObject, 'Value') == false
    
    handles.QCButton.Value = true;
    
end

pltSrcTrs

% --- Executes on button press in XCORButton.
function XCORButton_Callback(hObject, eventdata, handles)
% hObject    handle to XCORButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles = guihandles;

%this aligns all the traces within the fitting window
useVec    = logical(getappdata(gcf,'useVec'));
Traces    = getappdata(gcf,'Traces');
fw_start  = getappdata(gcf,'fw_start');
fw        = getappdata(gcf,'fitting_window');
lh        = getappdata(gcf,'lineHandles');
t         = getappdata(gcf,'t');
cwf0      = getappdata(gcf,'cwf');

Tclipped = Traces;

%clip the traces and zero the deleted ones
for k = 1:length(Traces)
    
    Tclipped(k).data = Traces(k).data(t > fw_start & t < fw(2));
        
end

handles.XCORButton.String = 'Running...';
drawnow

dt         = zeros(size(Traces));
dt(useVec) = dt_LSQR( Tclipped(useVec) );
dtsamples  = round(dt*Traces(1).sampleRate);

for k = 1:length(Traces)
    
    Traces(k).data = circshift(Traces(k).data, dtsamples(k));
    lh(k).YData    = Traces(k).data + k;

    setappdata(gcf,'cwf',k);
    
    add_to_source%remove the old
    add_to_source%add the shifted
    
    Traces(k).dt = dt(k);
    
end

setappdata(gcf,'cwf',cwf0);
setappdata(gcf, 'Traces', Traces);

pltSrcTrs

handles.XCORButton.String = 'XCOR';
drawnow


% --- Executes on button press in FlipWaveforms.
function FlipWaveforms_Callback(hObject, eventdata, handles)
% hObject    handle to FlipWaveforms (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles = guihandles;

%this aligns all the traces within the fitting window
Traces = getappdata(gcf,'Traces');
lh     = getappdata(gcf,'lineHandles');
cwf0   = getappdata(gcf,'cwf');

for k = 1:length(Traces)
    
    Traces(k).data = -1*Traces(k).data;
    lh(k).YData    = Traces(k).data + k;
    
    setappdata(gcf,'cwf',k);
    
    add_to_source%remove the old
    add_to_source%add the shifted
        
end

setappdata(gcf, 'cwf',cwf0);
setappdata(gcf, 'Traces', Traces);

pltSrcTrs


% --- Executes on button press in IntButton.
function IntButton_Callback(hObject, eventdata, handles)
% hObject    handle to IntButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles = guihandles;

%this aligns all the traces within the fitting window
Traces = getappdata(gcf, 'Traces');
lh     = getappdata(gcf, 'lineHandles');
cwf0   = getappdata(gcf, 'cwf');
t      = getappdata(gcf, 't');
stfw   = getappdata(gcf, 'fw_start');
fw     = getappdata(gcf, 'fitting_window');
level  = getappdata(gcf, 'level');

level = level + 1;
setappdata(gcf, 'level',level);

if level == 1
   
    handles.DiffButton.String = 'Velocity';
    handles.DiffButton.Enable = 'on';
    handles.IntButton.Enable  = 'off';
    
elseif level == 0
    
    handles.DiffButton.String = 'Acceleration';
    handles.IntButton.String  = 'Displacement';
    handles.DiffButton.Enable = 'on';
    handles.IntButton.Enable  = 'on';
    
end

for k = 1:length(Traces)
    
    Traces(k).data = cumsum(Traces(k).data);
    Traces(k).data = Traces(k).data - mean(Traces(k).data);
    Traces(k).data = Traces(k).data/rms(Traces(k).data(t > stfw & t < fw(2)));
    
    lh(k).YData    = Traces(k).data + k;
    
    setappdata(gcf,'cwf',k);
    setappdata(gcf, 'Traces', Traces);
    
    add_to_source%remove the old
    add_to_source%add the shifted
        
end

setappdata(gcf, 'cwf',cwf0);
setappdata(gcf, 'Traces', Traces);

pltSrcTrs

% --- Executes on button press in DiffButton.
function DiffButton_Callback(hObject, eventdata, handles)
% hObject    handle to DiffButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles = guihandles;

%this aligns all the traces within the fitting window
Traces = getappdata(gcf, 'Traces');
lh     = getappdata(gcf, 'lineHandles');
cwf0   = getappdata(gcf, 'cwf');
t      = getappdata(gcf, 't');
stfw   = getappdata(gcf, 'fw_start');
fw     = getappdata(gcf, 'fitting_window');
level  = getappdata(gcf, 'level');

level = level - 1;
setappdata(gcf, 'level',level);

if level == -1
   
    handles.DiffButton.Enable = 'off';
    handles.IntButton.String  = 'Velocity';
    handles.IntButton.Enable  = 'on';
    
elseif level == 0
    
    handles.DiffButton.String = 'Acceleration';
    handles.IntButton.String  = 'Displacement';
    handles.DiffButton.Enable = 'on';
    handles.IntButton.Enable  = 'on';
    
end

for k = 1:length(Traces)
    
    Traces(k).data = [diff(Traces(k).data); 0];
    Traces(k).data = Traces(k).data - mean(Traces(k).data);
    Traces(k).data = Traces(k).data/rms(Traces(k).data(t > stfw & t < fw(2)));

    lh(k).YData    = Traces(k).data + k;
    
    setappdata(gcf,'cwf',k);
    setappdata(gcf, 'Traces', Traces);
    
    add_to_source%remove the old
    add_to_source%add the shifted
        
end

setappdata(gcf, 'cwf',cwf0);
setappdata(gcf, 'Traces', Traces);

pltSrcTrs


%%%%%Code saved after button deleted
% --- Executes on button press in ShowHideButton.
function ShowHideButton_Callback(hObject, eventdata, handles)
% hObject    handle to ShowHideButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles = guihandles;

cwf      = getappdata(gcf, 'currentWf');
t        = getappdata(gcf, 't');
useVec   = logical(getappdata(gcf, 'useVec'));
plotmode = getappdata(gcf, 'plotmode');
midpoint = getappdata(gcf, 'midpoint');
xrange   = getappdata(gcf, 'xrange');
    
if strcmp(hOject.String, 'Hide Deleted')
    
    Traces   = getappdata(gcf, 'Traces');
    Traces_all = Traces;    
    Traces     = Traces(useVec);
    
    lh = getappdata(gcf, 'lineHandles');
    delete(lh)

    axes(handles.allWf_ax)
    hold on
    %plot all of the traces
    for k = 1:length(Traces)

        line_handles(k) = plot(t, Traces(k).data + k, 'k-', 'lineWidth', 1);

    end

    if plotmode == 1
    
        line_handles(cwf).LineWidth = 4;
        line_handles(cwf).Color     = 'b';
        
    elseif plotmode == 2
       
        line_handles(cwf).LineWidth = 4;
        line_handles(cwf).Color     = 'k';
        
    end
    
    xlabel('Time, s')
    xlim([(midpoint - xrange) (midpoint + xrange)])
    ylim([-1 (length(line_handles)+1)])

    setappdata(gcf,'lineHandles',line_handles);
    setappdata(gcf, 'Traces', Traces);
    setappdata(gcf, 'Traces_all', Traces_all);
    
    hObject.String = 'Show Deleted';
    
elseif strcmp(hOject.String, 'Show Deleted')
    
    Traces = getappdata(gcf, 'Traces_all');
    fw     = getappdata(gcf, 'fitting_window');
    
    Traces_all = Traces;    
    Traces     = Traces(useVec);
    
    lh = getappdata(gcf, 'lineHandles');
    delete(lh)

    axes(handles.allWf_ax)
    hold on
    %plot all of the traces
    for k = 1:length(Traces)

        line_handles(k) = plot(t, Traces(k).data + k, 'k-', 'lineWidth', 1);

    end

    if plotmode == 1
    
        line_handles(cwf).LineWidth = 4;
        line_handles(cwf).Color     = 'b';
        
    elseif plotmode == 2
       
        line_handles(cwf).LineWidth = 4;
        line_handles(cwf).Color     = 'k';
        
    end
    
    xlabel('Time, s')
    xlim([(midpoint - xrange) (midpoint + xrange)])
    ylim([-1 (length(line_handles)+1)])

    setappdata(gcf, 'lineHandles',line_handles);
    setappdata(gcf, 'Traces', Traces);
    setappdata(gcf, 'Traces_all', Traces_all);
    
    hObject.String = 'Hide Deleted';
    
end


% --- Executes on button press in WindowButton.
function WindowButton_Callback(hObject, eventdata, handles)
% hObject    handle to WindowButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



stfwh          = getappdata(gcf, 'stfwh');
stfwh_all      = getappdata(gcf, 'stfwh_all');
fwh            = getappdata(gcf, 'fwh');
fwh_all        = getappdata(gcf, 'fwh_all');
t              = getappdata(gcf, 't');

delete(stfwh); delete(stfwh_all);
delete(fwh);   delete(fwh_all);

axes(handles.currWf_ax)

handles.WindowButton.String = 'Pick start';
[ fw_start, ~]       = ginput(1);
handles.WindowButton.String = 'Pick ends';
[ fitting_window, ~] = ginput(2);
fitting_window       = sort(fitting_window);
handles.WindowButton.String = 'Move Fit Window';

fw_start          = max([ 0 fw_start]);
fitting_window(1) = max([ 0 fitting_window(1)]);
fitting_window(2) = min([ fitting_window(2) max(t) ]);

axes(handles.currWf_ax)
ylimits = ylim;
stfwh = plot( [ fw_start fw_start ], [ ylimits(1) ylimits(2) ], 'k--', 'LineWidth', 1.5);
setappdata(gcf, 'stfwh', stfwh);

fwh(1) = plot( [ fitting_window(1) fitting_window(1) ], [ ylimits(1) ylimits(2) ], 'k--', 'LineWidth', 0.5);
fwh(2) = plot( [ fitting_window(2) fitting_window(2) ], [ ylimits(1) ylimits(2) ], 'k--', 'LineWidth', 0.5);
setappdata(gcf, 'fwh', fwh);

axes(handles.allWf_ax)
ylimits = ylim;
stfwh_all = plot( [ fw_start fw_start ], [ ylimits(1) ylimits(2) ], 'k--', 'LineWidth', 1.5);
setappdata(gcf, 'stfwh_all', stfwh_all);
fwh_all(1) = plot( [ fitting_window(1) fitting_window(1) ], [ ylimits(1) ylimits(2) ], 'k--', 'LineWidth', 0.5);
fwh_all(2) = plot( [ fitting_window(2) fitting_window(2) ], [ ylimits(1) ylimits(2) ], 'k--', 'LineWidth', 0.5);
setappdata(gcf, 'fwh_all', fwh_all);

setappdata(gcf, 'fw_start', fw_start);
setappdata(gcf, 'fitting_window', fitting_window);

pltSrcTrs
   
function LowHz_Callback(hObject, eventdata, handles)
% hObject    handle to LowHz (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of LowHz as text
%        str2double(get(hObject,'String')) returns contents of LowHz as a double


% --- Executes during object creation, after setting all properties.
function LowHz_CreateFcn(hObject, eventdata, handles)
% hObject    handle to LowHz (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function HighHz_Callback(hObject, eventdata, handles)
% hObject    handle to HighHz (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of HighHz as text
%        str2double(get(hObject,'String')) returns contents of HighHz as a double


% --- Executes during object creation, after setting all properties.
function HighHz_CreateFcn(hObject, eventdata, handles)
% hObject    handle to HighHz (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function MinDeg_Callback(hObject, eventdata, handles)
% hObject    handle to MinDeg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of MinDeg as text
%        str2double(get(hObject,'String')) returns contents of MinDeg as a double


% --- Executes during object creation, after setting all properties.
function MinDeg_CreateFcn(hObject, eventdata, handles)
% hObject    handle to MinDeg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function MaxDeg_Callback(hObject, eventdata, handles)
% hObject    handle to MaxDeg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of MaxDeg as text
%        str2double(get(hObject,'String')) returns contents of MaxDeg as a double


% --- Executes during object creation, after setting all properties.
function MaxDeg_CreateFcn(hObject, eventdata, handles)
% hObject    handle to MaxDeg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on slider movement.
function FWslide_Callback(hObject, eventdata, handles)
% hObject    handle to FWslide (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

pltSrcTrs
plotMap

% --- Executes during object creation, after setting all properties.
function FWslide_CreateFcn(hObject, eventdata, handles)
% hObject    handle to FWslide (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes during object creation, after setting all properties.
function map_ax_CreateFcn(hObject, eventdata, handles)
% hObject    handle to map_ax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate map_ax

function map_ax_ButtonDownFcn(hObject, eventdata, handles)

if eventdata.Button == 1
   
    Traces   = getappdata(gcf, 'Traces');
    lh       = getappdata(gcf, 'lineHandles');
    cwf      = getappdata(gcf, 'currentWf');
    plotmode = getappdata(gcf, 'plotmode');
    inSrcVec = getappdata(gcf, 'inSrcVec');

    lh(cwf).Color     = 'k';
    
    if ~inSrcVec(cwf)
    
        lh(cwf).LineWidth = 1;
        
    end
    
    min_dist = 1e10;
    
    for k = 1:length(Traces)
       
        testdist = distance(Traces(k).latitude, Traces(k).longitude, ...
            eventdata.IntersectionPoint(2), eventdata.IntersectionPoint(1));
        
        if min_dist > testdist
           
            min_dist = testdist; 
            cwf = k;
            
        end
        
    end
    
    if plotmode == 1
    
        lh(cwf).LineWidth = 4;
        lh(cwf).Color     = 'b';
        
    elseif plotmode == 2
       
        lh(cwf).LineWidth = 4;
        lh(cwf).Color     = 'k';
        
    end
    
    setappdata(gcf, 'lineHandles', lh);
    setappdata(gcf, 'currentWf',   cwf);
    pltSrcTrs
    plotMap

end


% --- Executes on button press in RemoveButton.
function RemoveButton_Callback(hObject, eventdata, handles)
% hObject    handle to RemoveButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = guihandles;

t        = getappdata(gcf, 't');
useVec   = logical(getappdata(gcf, 'useVec'));
midpoint = getappdata(gcf, 'midpoint');
xrange   = getappdata(gcf, 'xrange');
inSrcVec = getappdata(gcf, 'inSrcVec');
Traces   = getappdata(gcf, 'Traces');

Traces   = Traces(useVec);

lh = getappdata(gcf, 'lineHandles');
delete(lh)

axes(handles.allWf_ax)
hold on
%plot all of the traces
for k = 1:length(Traces)
    
    if inSrcVec(k)
    
        line_handles(k) = plot(t, Traces(k).data + k, 'k-', 'LineWidth', 4);
        
    else
       
        line_handles(k) = plot(t, Traces(k).data + k, 'k-', 'LineWidth', 1);
        
    end
    
    line_handles(k).ButtonDownFcn = @ClickCWF;

end

line_handles(1).LineWidth = 4;
line_handles(1).Color     = 'b';

xlabel('Time, s')
xlim([(midpoint - xrange) (midpoint + xrange)])
ylim([-1 (length(line_handles)+1)])

setappdata(gcf, 'lineHandles',line_handles);
setappdata(gcf, 'Traces', Traces);
setappdata(gcf, 'currentWf', 1);
setappdata(gcf, 'useVec', ones(size(Traces)));
    
function ClickCWF(src,evnt)

lh       = getappdata(gcf, 'lineHandles');
cwf      = getappdata(gcf, 'currentWf');
plotmode = getappdata(gcf, 'plotmode');
inSrcVec = getappdata(gcf, 'inSrcVec');

if inSrcVec(cwf)

    lh(cwf).LineWidth = 4;
    
else
   
    lh(cwf).LineWidth = 1;
    
end
    
lh(cwf).Color     = 'k';

cwf = round(mean(src.YData));%traces are zero mean, but need to be rounded

if plotmode == 1
    
    lh(cwf).LineWidth = 4;
    lh(cwf).Color     = 'b';
    
elseif plotmode == 2
    
    lh(cwf).LineWidth = 4;
    lh(cwf).Color     = 'k';
    
end
    
setappdata(gcf, 'lineHandles', lh);
setappdata(gcf, 'currentWf',   cwf);

pltSrcTrs
plotMap

function ClickMap(src,eventdata)

if eventdata.Button == 1
   
    Traces   = getappdata(gcf, 'Traces');
    lh       = getappdata(gcf, 'lineHandles');
    cwf      = getappdata(gcf, 'currentWf');
    plotmode = getappdata(gcf, 'plotmode');
    inSrcVec = getappdata(gcf, 'inSrcVec');

    if inSrcVec(cwf)

        lh(cwf).LineWidth = 4;

    else

        lh(cwf).LineWidth = 1;

    end
    
    lh(cwf).Color     = 'k';
    
    min_dist = 1e10;
    
    for k = 1:length(Traces)
       
        testdist = distance(Traces(k).latitude, Traces(k).longitude, ...
            eventdata.IntersectionPoint(2), eventdata.IntersectionPoint(1));
        
        if min_dist > testdist
           
            min_dist = testdist; 
            cwf = k;
            
        end
        
    end
    
    if plotmode == 1
    
        lh(cwf).LineWidth = 4;
        lh(cwf).Color     = 'b';
        
    elseif plotmode == 2
       
        lh(cwf).LineWidth = 4;
        lh(cwf).Color     = 'k';
        
    end
    
    setappdata(gcf, 'lineHandles', lh);
    setappdata(gcf, 'currentWf',   cwf);
    pltSrcTrs
    plotMap

end

function allWf_ax_ButtonDownFcn(hObject, eventdata, handles)

lh       = getappdata(gcf, 'lineHandles');

if round(eventdata.IntersectionPoint(2)) <= length(lh) && round(eventdata.IntersectionPoint(2)) >= 1

    cwf      = getappdata(gcf, 'currentWf');
    plotmode = getappdata(gcf, 'plotmode');
    inSrcVec = getappdata(gcf, 'inSrcVec');

    if inSrcVec(cwf) && plotmode == 1

        lh(cwf).LineWidth = 4;

    else

        lh(cwf).LineWidth = 1;

    end

    lh(cwf).Color     = 'k';

    cwf = round(eventdata.IntersectionPoint(2));%tracesd are zero mean, but need to be rounded

    if plotmode == 1

        lh(cwf).LineWidth = 4;
        lh(cwf).Color     = 'b';

    elseif plotmode == 2

        lh(cwf).LineWidth = 4;
        lh(cwf).Color     = 'k';

    end

    setappdata(gcf, 'lineHandles', lh);
    setappdata(gcf, 'currentWf',   cwf);

    pltSrcTrs
    plotMap

end
