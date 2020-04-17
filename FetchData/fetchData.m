function [Traces] = fetchData(E, S, pre, post, phase, tag, taper_fraction, sample_rate, ...
    channel_string, email, password)

myTrace = [];

% now double loop through events and stations

%seconds in fractional days (because that's what MATLAB uses for date/times)
pre  = pre/(24*60*60); 
post = post/(24*60*60);

for ke=1:length(E)
    
    disp([ num2str(ke) ' of ' num2str(length(E)) ]);
    
    Traces=[];
    for ks=1:length(S)
        
        %first check to see if this station was active when the earthquake
        %happened
        
        sDate=datenum(S(ks).StartDate);
        if isempty(S(ks).EndDate); S(ks).EndDate='2500-01-01 00:00:00.000'; end
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
        
        if phase == 'P'
        
            phase_time = pTime;
        
        elseif phase == 'S'
            
            phase_time = sTime;
            
        end
        
        if isempty(phase_time); phase_time=times(1); end
        
        originTimeStr  = E(ke).PreferredTime;
        originTimeNum  = datenum(originTimeStr);
        
        startTime = originTimeNum+phase_time/(24*60*60)-pre;
        startTime = datestr(startTime,'yyyy-mm-dd HH:MM:SS.FFF');
        
        endTime   = originTimeNum+phase_time/(24*60*60)+post;
        endTime   = datestr(endTime,'yyyy-mm-dd HH:MM:SS.FFF');
        
        %fprintf('Trying event #%.0f station %s\n',ke,S(ks).StationCode)
        try
            
            if strcmp(email, '')
            
                myTrace = irisFetch.Traces(S(ks).NetworkCode,...
                        S(ks).StationCode,'*',channel_string,startTime, endTime);
                    
            else
               
                myTrace = irisFetch.Traces(S(ks).NetworkCode,...
                        S(ks).StationCode,'*',channel_string,startTime, endTime, ...
                        { email, password });
                
            end
                
        catch ME

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
            
            if phase == 'P'
            
                [myTrace.phaseTimes]=deal([0 (sTime-pTime)]+pre*(24*60*60)); %pre back to seconds
            
            elseif phase == 'S'
               
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
    data = {Traces.data};
    tf_empty  = cellfun('isempty',data);
    Traces    = Traces(~tf_empty);
    
    for q = 1:length(Traces)
    
        %demean/detrend
        Traces(q).data = detrend(Traces(q).data);
        %taper
        Traces(q).data = tukeywin(length(Traces(q).data), taper_fraction).*Traces(q).data;
    
    end
        
    Traces   = wfResample(Traces, sample_rate);
    Traces   = wfRemInstResp(Traces);
    
    eventData = E(ke);
    fname     = sprintf('%s_%03.0f',tag,ke);
    
    save(fname,'Traces','eventData')
    
end