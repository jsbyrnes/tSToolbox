%resampleData

clear

block='./TA/';
%orid={ '042'};

files = dir([ block 'T*.mat']);

load('./MetaData/TARio');

for k = 1:length(files)
    
    fname1=[block files(k).name];
    
    load(fname1)
           
    %resampling everything to 40 sps    
    Traces = wfResample(Traces,40);
    Traces_pre = Traces;
    
    Traces = wfRemInstResp(Traces);
    
    name_tmp = files(k).name;
    
    name_tmp = [ name_tmp(1:end-4) 'corr.mat' ];
        
    %save in new mat file
    eventData = E(k);
    save( [ block name_tmp ],'Traces','eventData')
    
    fprintf('done with %s\n',[ block name_tmp ])
    
end