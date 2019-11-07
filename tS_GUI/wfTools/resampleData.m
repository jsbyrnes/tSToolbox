%resampleData

clear

block='./';
%orid={ '042'};

samplerate = 20;
files      = dir([ block 'TARio_118*.mat']);

for k = 1:length(files)
    
    fname1=[block files(k).name];
    
    load(fname1)
    
    for q = 1:length(Traces)
    
        %demean/detrend
        Traces(q).data = detrend(Traces(q).data);
        %taper
        Traces(q).data = tukeywin(length(Traces(q).data), 0.8).*Traces(q).data;
    
    end
        
    %resampling everything to 40 sps    
    Traces   = wfResample(Traces, samplerate);
    Traces   = wfRemInstResp(Traces);
    
    name_tmp = files(k).name;
    
    name_tmp = [ name_tmp(1:end-4) 'corr.mat' ];
        
    %save in new mat file
    %eventData = E(k);
    save( [ block name_tmp ],'Traces')
    
    fprintf('done with %s\n',[ block name_tmp ])
    
end