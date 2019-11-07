function [S1] = estimateSource_tmp(traces)

%Normalize the traces

for k=1:length(traces)
    
    traces(k).data=traces(k).data/max(traces(k).data);
        
end

%Align and Stack

%%%%%%%hardwiring here
% pre  = 10*traces(1).sampleRate;
% post = 10*traces(1).sampleRate;
% 
% AlignedTraces=NaN(pre+post+1,length(traces));
% 
for k=1:length(traces)
   
    %[m ix]=max(traces(k).data(1:150)); %<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< min or max
        %find min or max for each trace. Remove (1:150)? 
    %AlignedTraces(:,k)=traces(k).data(ix-pre:ix+post);
    %JSB edit aligned going in but need this to stay so I can keep the
    %fitting window fixed later in the code    
    AlignedTraces(:,k)=traces(k).data;
      %save the data for each trace with set pre/post around max/min 
  
end

Stack=sum(AlignedTraces,2)/length(traces);
%find average for all aligned traces

%now integrate, taper, and differential to ensure it sums to zero
% Stack = cumtrapz(S);
% Stack = tukeywin(length(S), 0.95);
% Stack = [0 diff(Stack)];
% 
waveform=Stack;
waveform=waveform/max(abs(waveform));
%normalize waveform wrt max/min

figure(1); clf; hold on

for k=1:size(AlignedTraces,2)
       
    plot(AlignedTraces(:,k),'-','Color',[.4 .4 .4])
           
end

plot(waveform,'k-','LineWidth',2)

x1=1; x2=size(AlignedTraces,1);

xlim([x1 x2])
ylim([-1.1 1.1])

KeepPicking=true;

PointSet=[];
h1=[];h2=[];

hwx=20; hwy=.1;
delbx=x2-50;
delby=.8;

donbx=x1+50;
donby=.8;

DeleteBoxX=[delbx-hwx delbx+hwx delbx+hwx delbx-hwx delbx-hwx];
DeleteBoxY=[delby-hwy delby-hwy delby+hwy delby+hwy delby-hwy];

DoneBoxX=[donbx-hwx donbx+hwx donbx+hwx donbx-hwx donbx-hwx];
DoneBoxY=[donby-hwy donby-hwy donby+hwy donby+hwy donby-hwy];


plot(DeleteBoxX,DeleteBoxY,'r-','LineWidth',6)
plot(DoneBoxX,DoneBoxY,'b-','LineWidth',6)

while KeepPicking
    
    point=ginput(1);
    
    if inpolygon(point(1),point(2),DoneBoxX,DoneBoxY);
        KeepPicking=false;
        break
        
    elseif inpolygon(point(1),point(2),DeleteBoxX,DeleteBoxY);
                
        PointSet(end,:)=[];
        
    else
        
        
        pointDiff=[point(1)-(1:length(waveform))' (point(2)-waveform)*400];
        
        [m ix]=min(hypot(pointDiff(:,1),pointDiff(:,2)));
        
        if m < 6
            point=[ix waveform(ix)];
        end
        
        if abs(point(2)) < 0.025; point(2)=0; end
        
        PointSet=[PointSet; point];
        
    end
    
    PointSet=unique(PointSet,'rows','stable');
    
    PointSet0=[1 0; PointSet; length(waveform) 0];
    

    S=interp1(PointSet0(:,1),PointSet0(:,2),1:length(waveform),'pchip');

    
    if ishandle(h1); delete(h1); end
    if ishandle(h2); delete(h2); end
    
    h1=plot(PointSet(:,1),PointSet(:,2),'ro','MarkerFace','r','MarkerSize',6);
    h2=plot(S,'b-','LineWidth',2);
    
    disp(sum(S))
end

%%
Displacement=cumtrapz(S);

ix=find(abs(Displacement)>.1,1,'last');

Taper = ones(size(Displacement));

Taper(ix-10:ix+10)=linspace(1,0,21).^.2;

Taper(ix+10:end)=0;

try
    Displacement=Displacement.*Taper;
catch
    keyboard
end

S1=diff(Displacement);

figure(2); clf; hold on;
subplot(2,1,1); plot(S1)
subplot(2,1,2); plot(Displacement)






