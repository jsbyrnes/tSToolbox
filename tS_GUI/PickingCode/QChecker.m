function [] = QChecker(fname)
% [] = QChecker(block,orid)
% block, orid, are strings
% HELP IS OBSOLETE
% QChecker: Script for quality checking of synthetics.
% to use: 
% make sure the appropriate file names are in lines 19,20 and 72.
% run the script; you will see the observed and synthtetic waveforms for
% the first trace and two boxes (red and green) in the upper corners. There
% will be a '*' inside one of the boxes. The green box means that record is
% OK, the red box means it's not OK and should be taken out. 
% All events are initially OK. to change to the red box simply click inside
% the red box, to change back to the green click inside it. 
% When satisfied with the status of a trace click outside the boxes (but 
% inside the plot) to advance to the next trace.
% the trace number and the total number of Traces are displayed in the
% lower left.
% you have to go through all the Traces in order for the program to save
% the changes.
% While you're looking at these, make a note of the trace numbers of the
% ones that need to be flipped (low tech solution, but will work).

%fname=['./results/' block orid];
%fname = [ './CombinedResults/MAGIC_P_' orid 'rs_3Hz_result' ];

load([ fname '.mat' ])

% create a QC field on the Traces structure; all true
tt=true(size(Traces));
ttc=num2cell(tt);
[Traces(:).QC]=ttc{:};

inPy=[.6 .6 1 1 .6];

outPy=[.6 .6 1 1 .6];

for k=1:length(Traces)
    figure(1); clf; hold on
    
    plot(Traces(k).data/max(Traces(k).data),'b-','LineWidth',2)
    plot(ts_run(k).data/max(ts_run(k).data),'r-','LineWidth',2)
    
    ylim([-1 1.1])
    
    traceLength=length(Traces(k).data);
    
    
    rectHWidth=traceLength*.1; %rectangle half-width
    inCenterX=traceLength-2*rectHWidth;
    outCenterX=1+2*rectHWidth;
    
    inPx=inCenterX+rectHWidth*[-1 1 1 -1 -1];
    outPx=outCenterX+rectHWidth*[-1 1 1 -1 -1];
    
    plot(inPx,inPy,'g-','LineWidth',4)
    plot(outPx,outPy,'r-','LineWidth',4)
    
    %estr=[orid]; % event string
    sstr=Traces(k).station; %station tstring
    tStarStr = num2str(ts_run(k).tStar_WF);%num2str(colorby.tStar(k));
    kount=sprintf('%d/%d',k,length(Traces));
    %text(outCenterX-rectHWidth,-.8,[sstr ' ' tStarStr ' '  kount ' ' estr ],'fontsize',14,'fontweight','bold')
    text(outCenterX-rectHWidth,-.8,[sstr ' ' tStarStr ' '  kount ],'fontsize',14,'fontweight','bold')
    
    changing=true;
    
    while changing
        
        if Traces(k).QC
            h=plot(inCenterX,.8,'g*','MarkerSize',14,'LineWidth',4);
        else
            h=plot(outCenterX,.8,'r*','MarkerSize',14,'LineWidth',4);
        end
        
        a=ginput(1);
        
        if inpolygon(a(1),a(2),inPx,inPy)
            Traces(k).QC=true; delete(h)
        elseif inpolygon(a(1),a(2),outPx,outPy)
            Traces(k).QC=false; delete(h)
        else
            changing=false;
        end
        
    end
end

fnameOut=[fname 'QC.mat'];
save(fnameOut,'Traces')


close(1)