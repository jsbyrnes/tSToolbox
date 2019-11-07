function [Traces] = fetchData(E,S,varargin)
%this function goes and fetches data from the IRIS DMC using irisFetch. 
%USAGE: [Traces] = fetchData(E,S)
%       [Traces] = fetchData(E,S,fnameRoot)
% where E and S are event and station structures as would be produced by
% irisFetch.
% The function fetches, for each event, all traces recorded by stations in
% the station list, if any.
% The data is saved in a .mat file named traces### where ### is sequential
% numbers.
% The user can optionally specify a file name root with an optional third
% argument (see usage above).
% This version does multiple traces for each station


%set fname root to default and override if needed.
tag='traces';
if length(varargin)>=1
    tag=varargin{1};
end

myTrace = [];

if length(varargin)>=2
   
    PorS = varargin{2};
    PorS = upper(PorS);
    
    if PorS ~= 'P' && PorS ~= 'S'
        
        PorS = 'P';
        
    end
    
else
    
    PorS = 'P';
    
end

% now double loop through events and stations

pre  = 600/(24*60*60); %seconds in fractional days (because that's what MATLAB uses for date/times)
post = 600/(24*60*60);%1 hour after to capture all phases

for ke=1:length(E)
    
    disp([ num2str(ke) ' of ' num2str(length(E)) ]);
    
    Traces=[];
    for ks=1:length(S)
        
        %first check to see if this station was active when the earthquake
        %happened
        
        sDate=datenum(S(ks).StartDate);
        eDate=datenum(S(ks).EndDate);
        eqDate=datenum(E(ke).PreferredTime);
        
        if ~(eqDate>sDate && eqDate<eDate) %if eq not between eDate and sDate
            %disp('skipped one')
            continue
        end
        
        %calculate Delta
        [D,Az]=distance(E(ke).PreferredLatitude,E(ke).PreferredLongitude,S(ks).Latitude,S(ks).Longitude);
        
        urlstr=sprintf('http://service.iris.edu/irisws/traveltime/1/query?distdeg=%f&evdepth=%f&noheader=true&mintimeonly=true',D,E(ke).PreferredDepth);

        urlstr2=sprintf('http://service.iris.edu/irisws/traveltime/1/query?distdeg=%f&evdepth=%f&mintimeonly=true',D,E(ke).PreferredDepth);
        
        tStr=urlread(urlstr);

        tCell=textscan(tStr,'%f %f %s %f %f %f %f %f %s %s');
        
        phases=tCell{3};
        times=tCell{4};
        
        pTime=times(strcmpi(phases,'P')); %this is time to P arrival.
        sTime=times(strcmpi(phases,'S'));
        
        pTime=min(pTime);
        sTime=min(sTime);
        
        if PorS == 'P'
        
            phase_time = pTime;
        
        elseif PorS == 'S'
            
            phase_time = sTime;
            
        end
        
        if isempty(phase_time); phase_time=times(1); end
        
        originTimeStr=E(ke).PreferredTime;
        originTimeNum=datenum(originTimeStr);
        
        startTime=originTimeNum+phase_time/(24*60*60)-pre;
        startTime=datestr(startTime,'yyyy-mm-dd HH:MM:SS.FFF');
        
        endTime=originTimeNum+phase_time/(24*60*60)+post;
        endTime=datestr(endTime,'yyyy-mm-dd HH:MM:SS.FFF');
        
        %fprintf('Trying event #%.0f station %s\n',ke,S(ks).StationCode)
        try
            %channelList='BH1,BH2,BHZ,BHN,BHE,HH1,HH2,HHZ,HHE,HHN';
            %channelList='BH1,BH2,BHZ,BHN,BHE';
            channelList='BHZ,HHZ';
                myTrace = irisFetch.Traces(S(ks).NetworkCode,...
                    S(ks).StationCode,'*',channelList,startTime, endTime);
        catch ME
            %keyboard
            continue
        end
         
        
        %add phase names and times and add to the whole
        if ~isempty(myTrace)
            
            try
            [myTrace.phaseNames]=deal({'P','S'});
            catch ME
                keyboard
            end
            
            %note these times are relative to the start of the trace:
            
            if PorS == 'P'
            
                [myTrace.phaseTimes]=deal([0 (sTime-pTime)]+pre*(24*60*60)); %pre back to seconds
            
            elseif PorS == 'S'
               
                [myTrace.phaseTimes]=deal([ (pTime - sTime) 0]+pre*(24*60*60)); %pre back to seconds
                
            end
                
            %I know, this is growing inside  a loop.
            Traces = [Traces myTrace];
            
            assignin('base','T',Traces)
        end
        
        myTrace = [];
        
    end
    
    %skip to next event if this event had no traces whatsover
    if isempty(Traces)
        continue
    end
    
    % Get rid of empty structures:
    try
    data={Traces.data};
    catch ME
        keyboard
    end
    
    tf_empty=cellfun('isempty',data);
    
    Traces=Traces(~tf_empty);
    
    eventData=E(ke);
    fname=sprintf('%s_%03.0f',tag,ke);
    
    save(fname,'Traces','eventData')
    
end