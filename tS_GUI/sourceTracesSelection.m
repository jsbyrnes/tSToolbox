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

% Edit the above text to modify the response to  help sourceTracesSelection

% Last Modified by GUIDE v2.5 17-Apr-2020 13:55:01

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

set(gcf,'renderer','opengl') %this helps a lot with ginput sluggishness.

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

%set default parameters
setappdata(gcf, 'HighCorner', '3');
setappdata(gcf, 'LowCorner', '0.02');
setappdata(gcf, 'ChannelsKeep', 'BHZ');
setappdata(gcf, 'DeltaUpperLimit', '90');
setappdata(gcf, 'DeltaLowerLimit', '30');
setappdata(gcf, 'xrange', '20');

% set the keypress_fcn for the figure
set(gcf,'WindowKeyPressFcn',@chngTr)
set(gcf,'WindowScrollWheelFcn',@scrollFcn) %MB April 15 2020

%for moving the traces laterally by a variable amount
setappdata(gcf, 'trShiftX', 1);

% make the figure resizeable
set(gcf,'Resize','on')

% and now make the actual graphic objects resizable
f=gcf; c=f.Children;
for k=1:length(c)
    c(k).Units='normalized';
end

%some cosmetic stuff
handles.tSFWplot.Box='on';
handles.map_ax.Box='on';
handles.currWf_ax.Box='on';
handles.allWf_ax.Box='on';

% make the figure resizeable
set(gcf,'Resize','on')

% and now make the actual graphic objects resizable
f=gcf; c=f.Children;
for k=1:length(c)
    c(k).Units='normalized';
end

%some cosmetic stuff
handles.tSFWplot.Box='on';
handles.map_ax.Box='on';
handles.currWf_ax.Box='on';
handles.allWf_ax.Box='on';

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
    xrange         = str2num(getappdata(gcf,'xrange'));
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
    
    %initialize waitbar
    hWB=waitbar(0,'Getting Started','Name','t* grid search progress');
    
    for k = 1:length(fw_ending)
        
        %print message to command window
        fprintf([ 'On window ' num2str(k) ' of ' num2str(length(fw_ending)) ' candidate windows\n' ]);
        
        %update waitbar
        waitbar(k/length(fw_ending),hWB,[ 'On window ' num2str(k) ' of ' num2str(length(fw_ending)) ' candidate windows' ])
        
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
    delete(hWB);
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

        line_handles(k).ButtonDownFcn = @ClickCWF;
        
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
    setappdata(gcf, 'xrange', num2str(xrange));

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

%for browsing for a file
if ~isempty(key.Modifier) && strcmp(key.Modifier, 'control') && strcmp(key.Key, 'd')
   
    [fname, dirname]           = uigetfile('*.mat','Please select data file');
    pathStr                    = [dirname fname];
    handles                    = guihandles;
    handles.DataLoadBox.String = pathStr;
   
       
end

dataLoaded = getappdata(gcf, 'dataLoaded');

%this is for moving the traces laterally by more than one sample (up to 9)
k1=key.Key(1);
if any(k1=='123456789'); setappdata(gcf,'trShiftX',str2double(k1)*2);  end

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
                
                %this moves the star without replotting the entire map
                hMark=getappdata(gcf,'hMark');
                Traces=getappdata(gcf,'Traces');
                set(hMark,'XData',Traces(cwf).longitude);
                set(hMark,'YData',Traces(cwf).latitude);
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

                %this moves the star without replotting the entire map
                hMark=getappdata(gcf,'hMark');
                Traces=getappdata(gcf,'Traces');
                set(hMark,'XData',Traces(cwf).longitude);
                set(hMark,'YData',Traces(cwf).latitude);
            end

        elseif strcmp(key.Key,'leftarrow')

            cwf=getappdata(gcf,'currentWf');
            lh=getappdata(gcf,'lineHandles');

            Traces = getappdata(gcf,'Traces');
            Traces(cwf).data = circshift(Traces(cwf).data, -1*getappdata(gcf,'trShiftX'));        
            lh(cwf).YData = Traces(cwf).data + cwf;
            setappdata(gcf, 'Traces', Traces);
                       
            add_to_source%remove the old
            add_to_source%add the shifted
            setappdata(gcf,'trShiftX',1)% reset to move one sample at a time            
            
            pltSrcTrs

        elseif strcmp(key.Key,'rightarrow')

            cwf=getappdata(gcf,'currentWf');
            lh=getappdata(gcf,'lineHandles');

            Traces = getappdata(gcf,'Traces');
            Traces(cwf).data = circshift(Traces(cwf).data, 1*getappdata(gcf,'trShiftX'));        
            lh(cwf).YData = Traces(cwf).data + cwf;
            setappdata(gcf, 'Traces', Traces);

            add_to_source%remove the old
            add_to_source%add the shifted
            setappdata(gcf,'trShiftX',1)% reset to move one sample at a time            
            
            pltSrcTrs
            
        elseif strcmp(key.Key,'backspace')

            useFunc%flips if in or out of useVec
            pltSrcTrs
            plotMap
            
        elseif strcmp(key.Key,'space')

            add_to_source%flips if in or out of source
            pltSrcTrs
            
        end

    elseif ~isempty(key.Modifier)

        if strcmp(key.Modifier{1}, 'shift')%adjust the position of the window

            axes(handles.currWf_ax)
            midpoint = getappdata(gcf, 'midpoint');
            xrange   = str2num(getappdata(gcf, 'xrange'));
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
            setappdata(gcf,'xrange',   num2str(xrange)); %these are the traces that will go into the source estimate
            
        end
        
    end
    
    
    
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
        
            try
            
            srcTr(:, k) = srcTr(:, k)...
                /rms(srcTr(t > stfw & t < fw(2), k));
            
            catch
                
                keyboard
                
            end
            
        end
        
        hSrTrM=plot(t, srcTr,'-','color',[.6 .6 .6]);
        set(hSrTrM,'Tag','srcTrMany');
        
        MsrcTr=mean(srcTr,2);
        hSrTr=plot(t, MsrcTr,'-','color','k','linewidth',1.5);
        set(hSrTr,'Tag','srcTr');
        
        delete(findobj('Tag','theTrace')); %this just gets rid of the previous one
        hTr=plot(t, Traces(cwf).data/rms(Traces(cwf).data(t > stfw & t < fw(2))),'b', 'LineStyle', ls, 'linewidth',1.5);
        set(hTr,'Tag','theTrace')
        
        ylimits = ylim;
        stfwh = plot( [ stfw stfw ], [ ylimits(1) ylimits(2) ], 'k--', 'LineWidth', 1.5);
        setappdata(gcf, 'stfwh', stfwh);

        fwh(1) = plot( [ fw(1) fw(1) ], [ ylimits(1) ylimits(2) ], 'k--', 'LineWidth', 0.5);
        fwh(2) = plot( [ fw(2) fw(2) ], [ ylimits(1) ylimits(2) ], 'k--', 'LineWidth', 0.5);
        setappdata(gcf, 'fwh', fwh);
        
        %call to legend goes at the end, so it does not create new entries
        %after the window limits are drawn.
        h = legend([hSrTrM(1),hSrTr,hTr],'All source traces', 'Current source estimate', 'Current data','Location','southwest');
        set(h, 'Tag', 'srcLegend');
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
        
        hSynth=plot(t, ts_run(fw_ind, cwf).data/rms(ts_run(fw_ind, cwf).data(t > stfw & t < fw(2))),'r', 'LineStyle', ls, 'linewidth', 1.5);
        set(hSynth,'Tag','synTr');
        
        hObs=plot(t, Traces(cwf).data/rms(Traces(cwf).data(t > stfw & t < fw(2))),'-','color','k', 'LineStyle', ls, 'linewidth',1.5);
        set(hObs,'Tag','srcTr');
        
    end
    
    
        
    fw = getappdata(gcf, 'fitting_window');
    %fwh(2).XData   = [ (fw(1) + (handles.FWslide.Value/Traces(1).sampleRate)) (fw(1) + (handles.FWslide.Value/Traces(1).sampleRate)) ];

    ylimits = ylim;
    stfwh = plot( [ stfw stfw ], [ ylimits(1) ylimits(2) ], 'k--', 'LineWidth', 1.5);
    setappdata(gcf, 'stfwh', stfwh);
    fwh(1) = plot( [ fw(1) fw(1) ], [ ylimits(1) ylimits(2) ], 'k--', 'LineWidth', 0.5);
    fwh(2) = plot( [ (fw(1) + (handles.FWslide.Value/Traces(1).sampleRate)) (fw(1) + (handles.FWslide.Value/Traces(1).sampleRate)) ], [ ylimits(1) ylimits(2) ], 'k--', 'LineWidth', 0.5);
    fwh(1).Visible = 'off';
    setappdata(gcf, 'fwh', fwh);
    
    %call to legend goes at the end, so it does not create new entries
    %after the window limits are drawn.
    h = legend([hObs, hSynth],'Observed', [ 'Synthetic: \Deltat* = ' num2str(ts_run(fw_ind, cwf).tStar_WF, 3)],'Location','southwest');
    set(h, 'Tag', 'srcLegend');
    
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
xrange   = str2num(getappdata(gcf, 'xrange'));

try
xlim([ (midpoint - xrange) (midpoint + xrange) ]);
catch
keyboard
end

%ylim([-1 1]);
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
    hMapQC  = scatter( [ts_run(fw_ind, useVec).longitude],  [ts_run(fw_ind, useVec).latitude],  50, [ts_run(fw_ind, useVec).tStar_WF],  'filled');
    hMapNQC = scatter( [ts_run(fw_ind, ~useVec).longitude], [ts_run(fw_ind, ~useVec).latitude], 50, [ts_run(fw_ind, ~useVec).tStar_WF], 'filled');
    %colormap(cm); hC = colorbar; hC.Location = 'westoutside'; hC.Label.String = '\Deltat*, s';
    %caxis([ min([ts_run(fw_ind, useVec).tStar_WF] - 0.1) max([ts_run(fw_ind, useVec).tStar_WF] + 0.1)]);
    
    colormap(cm); 

    % edited by Zhao, fix the problem of multi colorbars when press
    % backspace (when changing the file path)
    C = findall(gcf,'type','ColorBar');
    if isempty(C)
        hC = colorbar; hC.Location = 'westoutside'; hC.Label.String = '\Deltat*, s';
        caxis([ min([ts_run(fw_ind, useVec).tStar_WF] - 0.1) max([ts_run(fw_ind, useVec).tStar_WF] + 0.1)]);
    end
    
    xlabel(['Longitude, ' char(176)]);
    ylabel(['Latitude, ' char(176)]);
    
    hMapQC.MarkerEdgeColor  = 'k';
    hMapNQC.MarkerFaceColor = 'k';
    
    hMapQC.ButtonDownFcn  = @ClickMap;
    hMapNQC.ButtonDownFcn = @ClickMap;
    
    hMark(1) = plot(ts_run(1, cwf).longitude, ts_run(1, cwf).latitude, 'k+', 'MarkerSize', 10, 'linewidth',1);
    hMark(2) = plot(ts_run(1, cwf).longitude, ts_run(1, cwf).latitude, 'wx', 'MarkerSize', 10, 'linewidth',1);

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
name     = handles.DataLoadBox.String;
dl       = getappdata(gcf, 'dataLoaded');

HighCorner      = str2num(getappdata(gcf, 'HighCorner'));
LowCorner       = str2num(getappdata(gcf, 'LowCorner'));
ChannelsKeep    = getappdata(gcf, 'ChannelsKeep');
DeltaUpperLimit = str2num(getappdata(gcf, 'DeltaUpperLimit'));
DeltaLowerLimit = str2num(getappdata(gcf, 'DeltaLowerLimit'));
xrange          = str2num(getappdata(gcf, 'xrange'));

if ~dl

    s = strsplit(name, '.');
    if ~strcmp('mat', s(end))

        name = [ name '.mat' ];%add the extension

    end
    
    %get the data and plot it all

    if exist(name, 'file')

        filter_bounds                   = [ (LowCorner) (HighCorner) ];
        filter_bounds(filter_bounds==0) = NaN;
        degVec = [ (DeltaLowerLimit) (DeltaUpperLimit) ];
        
        load(name)
        
        handles.LoadDataButton.String = 'Loading...';
        drawnow
        
        if ~exist('fw', 'var')%something that's saved
        
            if isempty(Traces)

                error('No traces with the given channel(s)');

            end

            lon = [Traces.longitude];
            [~, staind] = sort(lon);
            Traces = Traces(staind);

            % remove instrument response only when it has not been removed
            % during fetching
            if (any([Traces.instrument]))

                Traces = wfRemInstResp(Traces); 

            end

            if ~isnan(filter_bounds)

                Traces = wfButterworth( Traces, filter_bounds);

            end

            for k = 1:length(Traces)

                Traces(k).data = Traces(k).data - mean(Traces(k).data);
                Traces(k).data = Traces(k).data/rms(Traces(k).data);

            end

            % equalize trace length first
            minlen = 1e100;
            for k = 1:length(Traces)

                minlen = min([minlen length(Traces(k).data)]);

            end

            for k = 1:length(Traces)

                Traces(k).data        = Traces(k).data(1:minlen);
                Traces(k).sampleCount = minlen;

            end

            xd             = minlen/Traces(1).sampleRate;
            midpoint       = (xd/2);
            xlimits        = [ (midpoint - xrange) ((midpoint) + xrange) ];
            fitting_window = [ (midpoint - xrange/4) ((midpoint) + xrange/4) ];
            fw_start       = (midpoint - xrange/3);
            
            % set "use" to yes for all of them unless they are out of range
            useVec = true(1, length(Traces));

            for k = 1:length(Traces)

                [arc, azi] = distance(eventData.PreferredLatitude, eventData.PreferredLongitude, ...
                    Traces(k).latitude, Traces(k).longitude);

                if arc < degVec(1) || arc > degVec(2)

                    useVec(k) = 0;

                end

                Traces(k).evtarc = arc;
                Traces(k).evtazi = azi;

            end

            if (strcmp(ChannelsKeep, 'R') || strcmp(ChannelsKeep, 'T'))
                                           
                %rotate the two horizontals into R and T, then keep the one
                %you want, rename and store
                
                sta = unique({Traces.station});
                
                for j = 1:length(sta)
                   
                    %get horizontals
                    ind = find(strcmp({Traces.station}, sta{j}) & ([Traces.dip] == 0));
                    
                    if isempty(ind)
                        
                        continue
                        
                    end
                    
                    hT = Traces(ind);
                    
                    %now figure out which one is +90 of the other
                    azimuths = [hT.azimuth];
                                                            
                    if (length(azimuths) > 2)
                        
                        for qqq = 1:length(ind)
                        
                            Traces(ind(qqq)).channel = { 'null' };
                            
                        end
                        
                        continue
                        
                    end
                    
                    if wrapTo360(azimuths(1) - azimuths(2)) < wrapTo360(azimuths(2) - azimuths(1))
                        
                        aind = [ 2 1 ];
                        
                    else
                        
                        aind = [ 1 2 ];
                        
                    end
                    
                    azimuths = azimuths(aind);
                                        
                    if (abs(wrapTo360(diff(azimuths) - 90)) > 2) %slight tolerence for weird entries
                        
                        for qqq = 1:length(ind)
                        
                            Traces(ind(qqq)).channel = { 'null' };
                            
                        end
                        
                        continue
                        
                    end
                                        
                    hT  = hT(aind);
                    ind = ind(aind);
                    
                    a  = wrapTo360(hT(1).evtazi - azimuths(1));
                    
                    rt = [ cosd(a) sind(a); -sind(a) cosd(a) ]*[hT(1).data'; hT(2).data'];
                    
                    Traces(ind(1)).data = rt(1, :)';
                    Traces(ind(2)).data = rt(2, :)';
                    Traces(ind(1)).channel = 'R';
                    Traces(ind(2)).channel = 'T';
                                        
                end
                
            end
            
            chan = {Traces.channel};
            Traces(~strcmp(chan, ChannelsKeep)) = [];
            useVec(~strcmp(chan, ChannelsKeep)) = [];
            
            if ~any(useVec)
                
                f = msgbox('No traces within range', 'Error','error');
                a=f.Children(3); a.Children.FontSize=14;
                error('No traces within range');
                
            end
        
            % set "inSrc" to no for all of them
            inSrcVec = false(1, length(Traces));
            
            setappdata(gcf,'forSrc',        {[] []}); %these are the traces that will go into the source estimate
            setappdata(gcf,'ts_run',        []); %clear any preexisting_run
            setappdata(gcf,'level',         0); %clear any preexisting_run. NOT SAVED IN EARLIEST VERSIONS
            
            t = (1:minlen)/Traces(1).sampleRate;
            
            handles.QCButton.Enable            = 'off';
            
        else %you are reloading an old pick
            
            minlen         = length(Traces(1).data);
            xrange         = (minlen/Traces(1).sampleRate)/2;
            midpoint       = xrange;
            xlimits        = [ (midpoint - xrange) ((midpoint) + xrange) ];
            fitting_window = fw;%fw_start already loaded in correct name

            % set "use" to yes for all of them unless they are out of range
            useVec   = [Traces.QC]';
            inSrcVec = [Traces.inSrc]';

            t = (1:minlen)/Traces(1).sampleRate;
            
            handles.FWslide.Max         = diff(round((fitting_window)*Traces(1).sampleRate));
            handles.FWslide.SliderStep  = [ 1/handles.FWslide.Max 3/handles.FWslide.Max ];
            handles.FWslide.Value       = fw_ind;
            
            srcind = find(inSrcVec);
            for k = 1:length(srcind)
               
                srcTr(:, k) = Traces(srcind(k)).data;
                
            end
            
            setappdata(gcf,'forSrc', {srcind srcTr}); %these are the traces that will go into the source estimate
            setappdata(gcf,'ts_run', ts_run_all); 
            
            if exist('level', 'var')
            
                setappdata(gcf,'level',  level); %clear any preexisting_run. NOT SAVED IN EARLIEST VERSIONS
            
            else
                
                setappdata(gcf,'level',  0); %Will be incorrect if not saved from earlier run
                
            end
                
            handles.QCButton.Enable            = 'on';
            
        end
                    
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
                
        axes(handles.allWf_ax)
        hold on
        %plot all of the traces
        
        for k = 1:length(Traces)

            Traces(k).data = Traces(k).data/max(abs(Traces(k).data...
                (t > fw_start & t < fitting_window(2))));

            %variable name is changed because when you delete above, you
            %get an empty array. When you assign to the empty array, the
            %handle is a double, not an object. Seems a bug in matlab?
            
            if inSrcVec(k)
                
                lw = 4;
                
            else
                
                lw = 1;
                
            end
            
            if useVec(k)
            
                lh(k) = plot(t, Traces(k).data + k, 'k-', 'lineWidth', lw);
                                
            else
                
                lh(k) = plot(t, Traces(k).data + k, 'k--', 'lineWidth', 1);
                
            end

            lh(k).ButtonDownFcn = @ClickCWF;
            
        end

        xlabel('Time, s')
        xlim(xlimits)
        ylim([-1 (length(lh)+1)])
        
        % plot the first trace in the currWf axis
        setappdata(gcf,'currentWf',     1);
        setappdata(gcf,'nTraces',       length(Traces));
        setappdata(gcf,'lineHandles',   lh);
        setappdata(gcf,'Traces',        Traces);
        setappdata(gcf,'useVec',        useVec);
        setappdata(gcf,'inSrcVec',      inSrcVec);
        setappdata(gcf,'midpoint',      midpoint); %these are the traces that will go into the source estimate
        setappdata(gcf,'xrange',        num2str(xrange)); %these are the traces that will go into the source estimate
        setappdata(gcf,'t',             t);
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
cwf0      = getappdata(gcf,'currentWf');

Tclipped = Traces;

%clip the traces and zero the deleted ones
for k = 1:length(Traces)
    
    Tclipped(k).data = Traces(k).data(t > fw_start & t < fw(2));
        
end

handles.XCORButton.String = '...';
drawnow

dt         = zeros(size(Traces));
dt(useVec) = dt_LSQR( Tclipped(useVec) );
dtsamples  = round(dt*Traces(1).sampleRate);

for k = 1:length(Traces)
    
    Traces(k).data = circshift(Traces(k).data, dtsamples(k));
    try
    lh(k).YData    = Traces(k).data + k;
    catch
    keyboard
    end
    
    setappdata(gcf,'currentWf',k);
    
    add_to_source%remove the old
    add_to_source%add the shifted
    
    Traces(k).dt = dt(k);
    
end

setappdata(gcf,'currentWf',cwf0);
setappdata(gcf, 'Traces', Traces);

pltSrcTrs

handles.XCORButton.String = 'XC';
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
cwf0   = getappdata(gcf,'currentWf');

Traces(cwf0).data = -1*Traces(cwf0).data;
lh(cwf0).YData    = Traces(cwf0).data + cwf0;

setappdata(gcf, 'Traces', Traces);

add_to_source%remove the old
add_to_source%add the flipped waveform
        
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
cwf0   = getappdata(gcf, 'currentWf');
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
    
    setappdata(gcf,'currentWf',k);
    setappdata(gcf, 'Traces', Traces);
    
    add_to_source%remove the old
    add_to_source%add the shifted
        
end

setappdata(gcf, 'currentWf',cwf0);
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
cwf0   = getappdata(gcf, 'currentWf');
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
    
    setappdata(gcf,'currentWf',k);
    setappdata(gcf, 'Traces', Traces);
    
    add_to_source%remove the old
    add_to_source%add the shifted
        
end

setappdata(gcf, 'currentWf',cwf0);
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
xrange   = str2num(getappdata(gcf, 'xrange'));
    
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
xrange   = str2num(getappdata(gcf, 'xrange'));
inSrcVec = getappdata(gcf, 'inSrcVec');
Traces   = getappdata(gcf, 'Traces');
ts_run   = getappdata(gcf, 'ts_run'); %NOTE: MB added this, it wasn't updating ts_run March27_2020

%ts_run may not exist yet. It gets created after fitting waveforms or when clicking on the map before fitting.
%the following takes care of that, same as in plotMap.
if isempty(ts_run)
    ts_run = getappdata(gcf, 'Traces');
    for k = 1:length(ts_run)
        ts_run(k).tStar_WF = 0;
    end
end



%actually remove the ones for which usevec is false;
Traces   = Traces(useVec);
ts_run   = ts_run(:,useVec); %the colon is because ts_run may have several rows. It likely does if fitWaveforms has run.

%remove the traces so we can re-plot
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
setappdata(gcf, 'ts_run', ts_run);
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


% --- Executes on button press in pushbutton13.
function pushbutton13_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton13 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

defAns = { num2str(getappdata(gcf, 'HighCorner')) num2str(getappdata(gcf, 'LowCorner')) ...
    num2str(getappdata(gcf, 'DeltaLowerLimit')) num2str(getappdata(gcf, 'DeltaUpperLimit')) ...
    getappdata(gcf, 'ChannelsKeep'), num2str(getappdata(gcf, 'xrange')) };

prmpt = { 'High frequency filter corner, Hz', 'Low frequency filter corner, Hz',...
    'Minimum Delta', 'Maximum Delta', 'Channel, comma seperated', 'Width of plotting window, s'};

sPar = inputdlg(prmpt, 'Set Data Parameters', 1, defAns);

if ~isempty(sPar)
    
    setappdata(gcf,'HighCorner',      sPar{1});
    setappdata(gcf,'LowCorner',       sPar{2});
    setappdata(gcf,'DeltaLowerLimit', sPar{3});
    setappdata(gcf,'DeltaUpperLimit', sPar{4});
    setappdata(gcf,'ChannelsKeep',    sPar{5});
    setappdata(gcf,'xrange',          sPar{6});
    
end

% --- Executes on key press with focus on pushbutton13 and none of its controls.
function pushbutton13_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to pushbutton13 (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.CONTROL.UICONTROL)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in LineUp_btn.
function LineUp_btn_Callback(hObject, eventdata, handles)
% hObject    handle to LineUp_btn (see GCBO)
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

handles.LineUp_btn.String = '...';
drawnow


%find the min value in each clipped trace
for k = 1:length(Tclipped)
    [~,minIx(k)]=min(Tclipped(k).data); 
end

%get the shift in samples by subtracting the mean from the ixs
dtsamples=-round(minIx-mean(minIx));
dt=dtsamples/Traces(1).sampleRate;

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

handles.LineUp_btn.String = 'LU';
drawnow

function scrollFcn(h_obj,evt)
handles=guihandles;
nTraces=getappdata(gcf,'nTraces');
yl=get(handles.allWf_ax,'ylim');
ydist=yl(2)-yl(1);
ydist=ydist*.1; %will move by a third of current
if evt.VerticalScrollCount < 0
    ydist=min(ydist,yl(1)); %to prevent from going off the top
    yl=yl-ydist;
elseif evt.VerticalScrollCount > 0
    ydist=min(ydist,(nTraces+1-yl(2))); %to prevent from going off the bottom
    yl=yl+ydist;
end
set(handles.allWf_ax,'ylim',yl)


% --- Executes on button press in pushbutton15.
function pushbutton15_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton15 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%keep only the traces that were not culled out
Traces         = getappdata(gcf,'Traces');
useVec         = find(getappdata(gcf,'useVec'));
inSrcVec       = find(getappdata(gcf, 'inSrcVec'));
fw_start       = getappdata(gcf,'fw_start');
fitting_window = getappdata(gcf,'fitting_window');
t              = getappdata(gcf,'t');
xrange         = getappdata(gcf,'xrange');
midpoint       = getappdata(gcf,'midpoint');

window_ind     = find( (t > midpoint - xrange) & (t < midpoint + xrange));
t              = t(window_ind);
fw_start       = fw_start - t(1);
fitting_window = fitting_window - t(1);
t              = t - t(1);

%cut the data and the source
for i = 1:length(Traces)
    
    Traces(i).data        = detrend(Traces(i).data(window_ind));
    Traces(i).data        = Traces(i).data/rms(Traces(i).data(t > fw_start & t < fitting_window(2)));
    Traces(i).sampleCount = length(Traces(i).data);
    
end

setappdata(gcf, 'Traces', Traces);
setappdata(gcf, 'fw_start', fw_start);
setappdata(gcf, 'fitting_window', fitting_window);
setappdata(gcf, 't', t);

setappdata(gcf, 'forSrc', { [] [] });
cwf0 = getappdata(gcf,'currentWf');
for k = 1:length(inSrcVec)
       
    setappdata(gcf, 'currentWf', useVec(k));

    add_to_source;%remove, then add back in
    
end
setappdata(gcf, 'currentWf', cwf0);

lh = getappdata(gcf, 'lineHandles');
delete(lh)

axes(handles.allWf_ax)
hold on
%plot all of the traces
for k = 1:length(Traces)
        
    line_handles(k)               = plot(t, Traces(k).data/max(abs(Traces(k).data)) ...
        + k, 'k-', 'lineWidth', 1);
    line_handles(k).ButtonDownFcn = @ClickCWF;

    if ~any(k==useVec)
       
        line_handles(k).LineStyle = '--';
        
    end
    
    if any(k==inSrcVec)
       
        line_handles(k).LineWidth = 4;
        
    end
    
    if k == cwf0
       
        line_handles(k).LineWidth = 4;
        line_handles(k).Color     = 'b';
        
    end
    
end

xlabel('Time, s')
xlim([1 max(t)])
ylim([-1 (length(line_handles)+1)])

setappdata(gcf,'lineHandles',line_handles);
setappdata(gcf,'t',t);

fwh            = getappdata(gcf, 'fwh');

delete(fwh);

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

stfwh_all      = getappdata(gcf, 'stfwh_all');
fwh_all        = getappdata(gcf, 'fwh_all');

delete(stfwh_all);
delete(fwh_all);

axes(handles.allWf_ax)
ylimits = ylim;
stfwh_all = plot( [ fw_start fw_start ], [ ylimits(1) ylimits(2) ], 'k--', 'LineWidth', 1.5);
setappdata(gcf, 'stfwh_all', stfwh_all);
fwh_all(1) = plot( [ fitting_window(1) fitting_window(1) ], [ ylimits(1) ylimits(2) ], 'k--', 'LineWidth', 0.5);
fwh_all(2) = plot( [ fitting_window(2) fitting_window(2) ], [ ylimits(1) ylimits(2) ], 'k--', 'LineWidth', 0.5);
setappdata(gcf, 'fwh_all', fwh_all);

pltSrcTrs

