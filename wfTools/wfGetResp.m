function [resp]=wfGetResp(trace)
% this function retrieves the response file for a trace using
% irisDetch.Resp. 
% USAGE: [resp]=wfGetResp(trace)
% the input "trace" is a trace structure (for a single trace) in the
% irisFetch format
% the output "resp" is a structure with the fields:
% Amp, numZeroes, numPoles, Zeros, Poles.
% note, maybe using "include pz" in the call to irisFetch.traces, would save
% me going through all this trouble.

trace.location;%edit JSB
% 
% respFile=irisFetch.Resp(trace.network,trace.station,trace.location,trace.channel,...
%     trace.startTime,trace.endTime);

respFile=irisFetch.Resp(trace.network,trace.station,'*',trace.channel,...
    trace.startTime,trace.endTime);


%% parse it
% resp is a long string of characters that I now have to parse.
% 10 is the decimal number equivalent to line feed in ascii
lineFeed=find(respFile==char(10));

%I'm going to try to break this up into lines
nLines=length(lineFeed)-1;
C=cell(nLines,1);
for k=1:nLines
    
    C{k}=respFile(lineFeed(k)+1:lineFeed(k+1)-1);
    
end


%each cell of C contains one line of the RESP file
%% Now find out how many poles and zeros there are and extract their values
% I'm assuming that the location in the file (i.e. the line number) where
% this information appear is always the same. That's why aix, nzix and npix are
% hard-wired and they represent the line numbers for these three lines:
% B053F07     A0 normalization factor:               +5.71400E+08
% B053F09     Number of zeroes:                      2
% B053F14     Number of poles:                       5
aix=20;
nzix=22;
npix=23;

normFactor=C{aix}; %this gets out the line
normFactor=str2double(normFactor(end-11:end)); %this is just the number

nZeroes=C{nzix}; %this gets out the line
nZeroes=str2double(nZeroes(end-1:end)); %this extracts just the number I need
%lines containing the zeros
zix=nzix+3+(1:nZeroes); %these lines contain the information on the zeros
Zeros=C(zix); %get them out of the cell array
%now I use cellfun to get the numeric values out of the text string. This
%assumes that the format is always the same, so the numbers denoting the
%position in the string are hard-wired.
Zeros=cellfun(@(s) str2double(s(19:30))+1i*str2double(s(33:44)),Zeros);

%This works exactly as above.
nPoles=C{npix};
nPoles=str2double(nPoles(end-1:end));
%lines containing the poles
pix=npix+4+nZeroes+(1:nPoles);
Poles=C(pix);
Poles=cellfun(@(s) str2double(s(19:30))+1i*str2double(s(33:44)),Poles);
%%
%a(62)

%JSB commented this out
%fprintf('%+9.5e\n',normFactor)
%fprintf('%+9.5e  %+9.5e\n',[real(Zeros) imag(Zeros)]')
%fprintf('%+9.5e  %+9.5e\n',[real(Poles) imag(Poles)]')

resp.Amp=normFactor;
resp.Numzeroes=nZeroes;
resp.Numpoles=nPoles;
resp.Zeros=Zeros;
resp.Poles=Poles;

