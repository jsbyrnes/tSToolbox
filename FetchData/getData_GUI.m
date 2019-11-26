function varargout = getData_GUI(varargin)
% GETDATA_GUI MATLAB code for getData_GUI.fig
%      GETDATA_GUI, by itself, creates a new GETDATA_GUI or raises the existing
%      singleton*.
%
%      H = GETDATA_GUI returns the handle to a new GETDATA_GUI or the handle to
%      the existing singleton*.
%
%      GETDATA_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GETDATA_GUI.M with the given input arguments.
%
%      GETDATA_GUI('Property','Value',...) creates a new GETDATA_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before getData_GUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to getData_GUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help getData_GUI

% Last Modified by GUIDE v2.5 21-Nov-2019 11:41:12

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @getData_GUI_OpeningFcn, ...
                   'gui_OutputFcn',  @getData_GUI_OutputFcn, ...
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
% End initialization code - DO NOT EDIT


% --- Executes just before getData_GUI is made visible.
function getData_GUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to getData_GUI (see VARARGIN)

% Choose default command line output for getData_GUI
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes getData_GUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);

addpath('../tS_GUI/wfTools/');

%activate the figure toolbars
set(gcf,'ToolBar','figure')

%plot coastline
subplot(handles.map_ax); hold on;
C=load('coast');
plot(C.long,C.lat,'k-'); daspect([1 1 1]); ylim([-90 75])



%for testing, put some initial values on the etxt boxes

set(handles.minLat_etxt,'string','0');
set(handles.maxLat_etxt,'string','14');
set(handles.minLon_etxt,'string','-72');
set(handles.maxLon_etxt,'string','-62');

set(handles.sTime_etxt,'string','2004-01-01 00:00:00');
set(handles.eTime_etxt,'string','2006-12-31 00:00:00');


%set initial values for Station Search Parameters
setappdata(gcf,'detl','');
setappdata(gcf,'net','');
setappdata(gcf,'sta','');
setappdata(gcf,'loc','');
setappdata(gcf,'chan','BHZ,HHZ');

%set initial values for Event Search Parameters
setappdata(gcf,'minZ',0);
setappdata(gcf,'maxZ',750);
setappdata(gcf,'minM',5.5);
setappdata(gcf,'maxM',10);

%set initial values for download parameters
setappdata(gcf, 'tBefore'       , 600);
setappdata(gcf, 'tAfter'        , 600);
setappdata(gcf, 'Taper_Fraction', 0.2);
setappdata(gcf, 'SampleRate'    , 20);
setappdata(gcf, 'DownloadPhase' , 'P');
setappdata(gcf, 'FileTag'       , 'DefaultName');
setappdata(gcf, 'ChannelString' , 'BHZ,HHZ');
setappdata(gcf, 'Email'         , '');
setappdata(gcf, 'Password'      , '');

% --- Outputs from this function are returned to the command line.
function varargout = getData_GUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function sTime_etxt_Callback(hObject, eventdata, handles)
% hObject    handle to sTime_etxt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of sTime_etxt as text
%        str2double(get(hObject,'String')) returns contents of sTime_etxt as a double


% --- Executes during object creation, after setting all properties.
function sTime_etxt_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sTime_etxt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function eTime_etxt_Callback(hObject, eventdata, handles)
% hObject    handle to eTime_etxt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of eTime_etxt as text
%        str2double(get(hObject,'String')) returns contents of eTime_etxt as a double


% --- Executes during object creation, after setting all properties.
function eTime_etxt_CreateFcn(hObject, eventdata, handles)
% hObject    handle to eTime_etxt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function maxLat_etxt_Callback(hObject, eventdata, handles)
% hObject    handle to maxLat_etxt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of maxLat_etxt as text
%        str2double(get(hObject,'String')) returns contents of maxLat_etxt as a double


% --- Executes during object creation, after setting all properties.
function maxLat_etxt_CreateFcn(hObject, eventdata, handles)
% hObject    handle to maxLat_etxt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function minLat_etxt_Callback(hObject, eventdata, handles)
% hObject    handle to minLat_etxt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of minLat_etxt as text
%        str2double(get(hObject,'String')) returns contents of minLat_etxt as a double


% --- Executes during object creation, after setting all properties.
function minLat_etxt_CreateFcn(hObject, eventdata, handles)
% hObject    handle to minLat_etxt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function minLon_etxt_Callback(hObject, eventdata, handles)
% hObject    handle to minLon_etxt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of minLon_etxt as text
%        str2double(get(hObject,'String')) returns contents of minLon_etxt as a double


% --- Executes during object creation, after setting all properties.
function minLon_etxt_CreateFcn(hObject, eventdata, handles)
% hObject    handle to minLon_etxt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function maxLon_etxt_Callback(hObject, eventdata, handles)
% hObject    handle to maxLon_etxt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of maxLon_etxt as text
%        str2double(get(hObject,'String')) returns contents of maxLon_etxt as a double


% --- Executes during object creation, after setting all properties.
function maxLon_etxt_CreateFcn(hObject, eventdata, handles)
% hObject    handle to maxLon_etxt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in getSta_btn.
function getSta_btn_Callback(hObject, eventdata, handles)
% hObject    handle to getSta_btn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% get latLonBox values
minLat=str2double(get(handles.minLat_etxt,'string'));
maxLat=str2double(get(handles.maxLat_etxt,'string'));
minLon=str2double(get(handles.minLon_etxt,'string'));
maxLon=str2double(get(handles.maxLon_etxt,'string'));

% get start and end times, as well as other search params
srchParams.sTime=get(handles.sTime_etxt,'string');
srchParams.eTime=get(handles.eTime_etxt,'string');
srchParams.latLonBox=[minLat maxLat minLon maxLon];

srchParams.detl=getappdata(gcf,'detl');
srchParams.net=getappdata(gcf,'net');
srchParams.sta=getappdata(gcf,'sta');
srchParams.loc=getappdata(gcf,'loc');
srchParams.chan=getappdata(gcf,'chan');

% call the function 
set(gcf,'pointer','watch'); drawnow

[S] = findStations(srchParams);
set(gcf,'pointer','arrow'); drawnow

% make the plot
subplot(handles.map_ax);
h1=plot([S.Longitude],[S.Latitude],'bv');
%add a box showing the searched-for area
h2=plot([minLon maxLon maxLon minLon minLon],...
    [minLat minLat maxLat maxLat minLat],'r-','LineWidth',2);
%store handles in appdata (to delete later if desired)
HH=getappdata(gcf,'someHandles');
setappdata(gcf,'someHandles',[HH h1 h2]);

%store the station structure in appdata
setappdata(gcf,'SStruct',S);

% --- Executes on button press in rbbox_btn.
function rbbox_btn_Callback(hObject, eventdata, handles)
% hObject    handle to rbbox_btn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles=guihandles(gcf);

%do the rubber band thing
k = waitforbuttonpress; 
point1 = get(handles.map_ax,'CurrentPoint');    % button down detected 
finalRect = rbbox;                 % return figure units 
point2 = get(handles.map_ax,'CurrentPoint');    % button up detected 
point1 = point1(1,1:2);            % extract x and y 
point2 = point2(1,1:2); 

% get the min and max lat and lon
minLat=min(point1(1,2),point2(1,2));
maxLat=max(point1(1,2),point2(1,2));
minLon=min(point1(1,1),point2(1,1));
maxLon=max(point1(1,1),point2(1,1));

%put the values onto the etxt fields where they go
set(handles.minLat_etxt,'string',num2str(minLat));
set(handles.maxLat_etxt,'string',num2str(maxLat));
set(handles.minLon_etxt,'string',num2str(minLon));
set(handles.maxLon_etxt,'string',num2str(maxLon));


% --- Executes on button press in sv2WS_btn.
function sv2WS_btn_Callback(hObject, eventdata, handles)
% hObject    handle to sv2WS_btn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%this simply puts the structure with the station info in the workspace
S=getappdata(gcf,'SStruct');
E=getappdata(gcf,'EStruct');
assignin('base','S',S);
assignin('base','E',E);


% --- Executes on button press in clr_btn.
function clr_btn_Callback(hObject, eventdata, handles)
% hObject    handle to clr_btn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
dH=getappdata(gcf,'someHandles');
delete(dH);
setappdata(gcf,'someHandles',[])

% --- Executes on button press in staParams_btn.
function staParams_btn_Callback(hObject, eventdata, handles)
% hObject    handle to staParams_btn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

detl=getappdata(gcf,'detl');
net=getappdata(gcf,'net');
sta=getappdata(gcf,'sta');
loc=getappdata(gcf,'loc');
chan=getappdata(gcf,'chan');

defAns={detl,net,sta,loc,chan};

prmpt={'Detail Level','Network Code','Station Name','Location Code','Channel(s)'};

sPar=inputdlg(prmpt,'Set Station Search Params',1,defAns);

if ~isempty(sPar)
    setappdata(gcf,'detl',sPar{1});
    setappdata(gcf,'net',sPar{2});
    setappdata(gcf,'sta',sPar{3});
    setappdata(gcf,'loc',sPar{4});
    setappdata(gcf,'chan',sPar{5});
end


% --- Executes on button press in gtEvtR_btn.
function gtEvtR_btn_Callback(hObject, eventdata, handles)
% hObject    handle to gtEvtR_btn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% get radial coord values
Lat=str2double(get(handles.lt_etxt,'string'));
Lon=str2double(get(handles.ln_etxt,'string'));
minRad=str2double(get(handles.minR_etxt,'string'));
maxRad=str2double(get(handles.maxR_etxt,'string'));

% get start and end times, as well as other search params
srchParams.sTime=get(handles.sTime_etxt,'string');
srchParams.eTime=get(handles.eTime_etxt,'string');
srchParams.radialCoords=[Lat Lon maxRad minRad];

srchParams.minZ=getappdata(gcf,'minZ');
srchParams.maxZ=getappdata(gcf,'maxZ');
srchParams.minM=getappdata(gcf,'minM');
srchParams.maxM=getappdata(gcf,'maxM');


% call the function 
set(gcf,'pointer','watch'); drawnow

[E] = findEvents(srchParams);
set(gcf,'pointer','arrow'); drawnow

% make the plot
subplot(handles.map_ax);
h1=plot([E.PreferredLongitude],[E.PreferredLatitude],'r.');
%add a line showing the searched-for area
ang=0:360;
[ly1,lx1]=reckon(Lat,Lon,minRad,ang);

%this is to avoid breaks in the line
if max(diff(lx1))>300
   [m, ix]= max(diff(lx1));
   lx1=[lx1(1:ix) NaN lx1(ix+1:end)];
   ly1=[ly1(1:ix) NaN ly1(ix+1:end)];
end

[ly2,lx2]=reckon(Lat,Lon,maxRad,ang);
%this is to avoid breaks in the line
if max(diff(lx2))>300
   [m, ix]= max(diff(lx2));
   lx2=[lx2(1:ix) NaN lx2(ix+1:end)];
   ly2=[ly2(1:ix) NaN ly2(ix+1:end)];
end

h2=plot(lx1,ly1,'g-');
h3=plot(lx2,ly2,'g-');

%store handles in appdata (to delete later if desired)
HH=getappdata(gcf,'someHandles');
setappdata(gcf,'someHandles',[HH h1 h2 h3]);

%store the station structure in appdata
setappdata(gcf,'EStruct',E);

function lt_etxt_Callback(hObject, eventdata, handles)
% hObject    handle to lt_etxt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of lt_etxt as text
%        str2double(get(hObject,'String')) returns contents of lt_etxt as a double


% --- Executes during object creation, after setting all properties.
function lt_etxt_CreateFcn(hObject, eventdata, handles)
% hObject    handle to lt_etxt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ln_etxt_Callback(hObject, eventdata, handles)
% hObject    handle to ln_etxt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ln_etxt as text
%        str2double(get(hObject,'String')) returns contents of ln_etxt as a double


% --- Executes during object creation, after setting all properties.
function ln_etxt_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ln_etxt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pick_btn.
function pick_btn_Callback(hObject, eventdata, handles)
% hObject    handle to pick_btn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
pnt=ginput(1);

set(handles.lt_etxt,'string',num2str(pnt(2)));
set(handles.ln_etxt,'string',num2str(pnt(1)));

% --- Executes on button press in evtParams_btn.
function evtParams_btn_Callback(hObject, eventdata, handles)
% hObject    handle to evtParams_btn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

minZ=getappdata(gcf,'minZ');
maxZ=getappdata(gcf,'maxZ');
minM=getappdata(gcf,'minM');
maxM=getappdata(gcf,'maxM');

defAns={num2str(minZ),num2str(maxZ),num2str(minM),num2str(maxM)};

prmpt={'Minimum Depth','Maximum Depth','Minimum Magnitude','Maximum Magnitude'};

ePar=inputdlg(prmpt,'Set Event Search Params',1,defAns);

if ~isempty(ePar)
    setappdata(gcf,'minZ',str2double(ePar{1}));
    setappdata(gcf,'maxZ',str2double(ePar{2}));
    setappdata(gcf,'minM',str2double(ePar{3}));
    setappdata(gcf,'maxM',str2double(ePar{4}));
end


function maxR_etxt_Callback(hObject, eventdata, handles)
% hObject    handle to maxR_etxt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of maxR_etxt as text
%        str2double(get(hObject,'String')) returns contents of maxR_etxt as a double


% --- Executes during object creation, after setting all properties.
function maxR_etxt_CreateFcn(hObject, eventdata, handles)
% hObject    handle to maxR_etxt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function minR_etxt_Callback(hObject, eventdata, handles)
% hObject    handle to minR_etxt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of minR_etxt as text
%        str2double(get(hObject,'String')) returns contents of minR_etxt as a double


% --- Executes during object creation, after setting all properties.
function minR_etxt_CreateFcn(hObject, eventdata, handles)
% hObject    handle to minR_etxt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over staParams_btn.
function staParams_btn_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to staParams_btn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton9.
function pushbutton9_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

tB = getappdata(gcf,'tBefore');
tA = getappdata(gcf,'tAfter');
DP = getappdata(gcf,'DownloadPhase');
FT = getappdata(gcf,'FileTag');
TF = getappdata(gcf,'Taper_Fraction');
SR = getappdata(gcf,'SampleRate');
CS = getappdata(gcf,'ChannelString');
EM = getappdata(gcf,'Email');
PW = getappdata(gcf,'Password');

defAns = { num2str(tB), num2str(tA), DP, FT, num2str(TF), num2str(SR), CS, EM, PW};
prmpt  = {'Time before (seconds)', 'Time after (seconds)','Phase, P or S', ...
    'File Tags', 'Taper Fraction', 'Sample Rate to resample to', 'Channels (comma seperated, no spaces',...
    'Email (for restricted data)', 'Password (for restricted data)'};
dPar   = inputdlg(prmpt,'Set Download Parameters', 1, defAns);

if ~isempty(dPar) 
    
    setappdata(gcf, 'tBefore'       , str2double(dPar{1}));
    setappdata(gcf, 'tAfter'        , str2double(dPar{2}));
    setappdata(gcf, 'DownloadPhase' , dPar{3});
    setappdata(gcf, 'FileTag'       , dPar{4});
    setappdata(gcf, 'Taper_Fraction', str2double(dPar{5}));
    setappdata(gcf, 'SampleRate'    , str2double(dPar{6}));
    setappdata(gcf, 'ChannelString' , dPar{7});
    setappdata(gcf, 'Email'         , dPar{8});
    setappdata(gcf, 'Password'      , dPar{9});
    
end

% --- Executes on button press in pushbutton10.
function pushbutton10_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

E  = getappdata(gcf,'EStruct');
S  = getappdata(gcf,'SStruct');
tB = getappdata(gcf,'tBefore');
tA = getappdata(gcf,'tAfter');
DP = getappdata(gcf,'DownloadPhase');
FT = getappdata(gcf,'FileTag');
TF = getappdata(gcf,'Taper_Fraction');
SR = getappdata(gcf,'SampleRate');
CS = getappdata(gcf,'ChannelString');
EM = getappdata(gcf,'Email');
PW = getappdata(gcf,'Password');

dummy = fetchData(E, S, tB, tA, DP, FT, TF, SR, CS, EM, PW);
