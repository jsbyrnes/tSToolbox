%fit a smooth surface to all t* observations through an inversion scheme.
%This inversion is quick and useful for many applications. 
clear, close all
%% collect all lats, lons and t* by loading one results file at a time.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%EDIT THESE LINES TO CONTROL THE INVERSION
smoothness          = 5;
station_smallness   = 5;
model_smallness     = 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%EDIT THIS LINE TO LOAD THE DATA
fnames = dir('../MM_raw/*Measurement.mat');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%
%Save results as 
name = 'LinearInv2';
%%%%%%

%%%%%%
%These lines control how things are plotted
v                   = -0.15:0.01:0.15;
clim                = [-0.15 0.15];
save_for_gmt        = 0;%not possible in this verison of the code
rotation            = 55;
collapse_y          = 0;
buffer              = 20;%in km on the outside

label_for_colorbar  = '\Deltat*, s';

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Minor inversion parameters
nodeSpacing         = 10;
maxDist             = nodeSpacing*5;
QC                  = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%LEAVE THIS LINE UNCHANGED UNLESS YOU KNOW WHY YOU ARE CHANGING IT
fieldtoload = 'tStar_WF';
%%%%%%%%%%%%%%%%%%%%%%%%%

%% Code below here is the actual inversion
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
allLats   = [];
allLons   = [];
allTS     = [];
allSig    = [];
dataE     = [];
allSta    = {};

if isempty(fnames)
    
    error('No files found')
    
end

for k=1:length(fnames)
    
    fname=[ fnames(k).folder '\' fnames(k).name ];
    
    load(fname, 'ts_run', 'Traces');
    
    if QC
       
        use = [Traces.QC];
        
        ts_run = ts_run(use);
        
    end
    
    %loop to remove doubled lats/lons (not my preferred solution)
    for kt=1:length(ts_run)
        
        lt=ts_run(kt).latitude;
        ts_run(kt).latitude=lt(1);
        
        ln=ts_run(kt).longitude;
        ts_run(kt).longitude=ln(1);
        
    end
    
    sta   = {ts_run.station};
    lats  = [ts_run.latitude];
    lons  = [ts_run.longitude];
    tS    = [ts_run.(fieldtoload)];
    sig   = 0.3*ones(size([ts_run.(fieldtoload)]));%[traces.sigma];
    
    allLats  = [ allLats lats ];
    allLons  = [ allLons lons ];
    allTS    = [ allTS tS ];
    allSig   = [ allSig sig ];
    allSta   = [ allSta sta ];
    dataE    = [ dataE lats*0+k ]; % this is just the event number
    
end

uSta=unique(allSta); %list of all station names

for k = 1:length(unique(dataE))
    
    allTS(dataE==k) = allTS(dataE==k) - mean(allTS(dataE==k));
    
end

%% define a map projection and inversion grid
%origin in the center of the array footprint

centerLat = min(allLats)+(max(allLats)-min(allLats))/2;
centerLon = min(allLons)+(max(allLons)-min(allLons))/2;

origin = [ centerLat centerLon ];           % array center [ lat lon ]

mstruct = defaultm('mercator');
mstruct.origin = [ origin rotation ];%second number is rotation
mstruct = defaultm( mstruct );
mstruct.scalefactor = 6371;

%%%%%%%%%%
%%%consider rotating the grid for the inversion so cut down the number
%%%of unsed nodes
%%%%%%%%%%

%projected lat/lon of observations to X/Y:
[dataX, dataY] = mfwdtran(mstruct,allLats,allLons);
%figure; plot(dataX,dataY,'k.');

if collapse_y
   
    dataY = zeros(size(dataY));
    
end

%make the grid on which to invert
minX=min(dataX)-buffer;
maxX=max(dataX)+buffer;
minY=min(dataY)-buffer;
maxY=max(dataY)+buffer;

xVec=minX:nodeSpacing:maxX;
yVec=minY:nodeSpacing:maxY;

[xMat, yMat]=meshgrid(xVec,yVec);

%hold on

%plot(xMat(:),yMat(:),'rs')


%% Let's invert
% zeroeth step: sizes of things
nData=length(allTS); %number of observations
nModel=numel(xMat);  %number of model nodes
nEvents=length(fnames); %number of events
nSta=length(uSta);

% and change variable names:
dataV=allTS';
dataI=dataV*0;

for k=1:length(dataI)
    
   dataI(k)=find(strcmp(allSta{k},uSta));

end

% First make the G, one row at a time:
% x and y of each model point as row vectors
modelX=xMat(:); modelX=modelX';
modelY=yMat(:); modelY=modelY';

%maximum distance away from model node for which observations are relevant

%this is the power to which we elevate the inverse of the distance for
%weighting.
dPow = 1;

%allocate G
G = nan(nData,nModel);

for k = 1:nData
    
    stationX = dataX(k);
    stationY = dataY(k);
    
    sta_m_dist=sqrt( (modelX-stationX).^2 + (modelY-stationY).^2 );

    Grow = 1./sta_m_dist;
    
    ix       = (sta_m_dist>maxDist);
    Grow(ix) = 0;
    
    Grow = Grow.^dPow;
    %now normalize the row of G
    Grow   = Grow/sum(Grow);   
    G(k,:) = Grow;
    
end

%event term
%allocate
Gevt = nan(nData,nEvents);
%build row-by-row
for k = 1:nData
   Gevt(k,:)        = 0;
   Gevt(k,dataE(k)) = 1;
end

%station term
%allocate
Gsta = nan(nData,nSta);
%build row-by-row
for k = 1:nData
   Gsta(k,:)        = 0;
   Gsta(k,dataI(k)) = 1;
end

%Put sensitivity matrix together:
G = [G Gevt Gsta];

%station smallness constraint
GSsml = eye(nSta);
GSsml = GSsml*station_smallness; % smallness control factor

%model smallness
GMsml = eye(nModel);
GMsml = GMsml*model_smallness; % smallness control factor

%add regularization matrices
Gexpanded = [G; [zeros(nSta,nModel) zeros(nSta,nEvents) GSsml]];
Gexpanded = [Gexpanded; [GMsml zeros(nModel,nEvents) zeros(nModel,nSta)]];

%smoothness constraint
makeSmoothness % this script produces Gsmth
Gsmth     = Gsmth*smoothness;
Gexpanded = [Gexpanded; [Gsmth zeros(nModel,nEvents) zeros(nModel,nSta)]];

%% build expanded d
dExpanded = [dataV; zeros(nSta,1); zeros(nModel,1); zeros(nModel,1)];
%the added zeros are for station term smallness, model smallness and model
%smoothness.

%% and invert
Minverted = Gexpanded\dExpanded;% lsqr(Gexpanded, dExpanded, 1e-10, 3000);%
MI        = Minverted;

MImodel             = Minverted(1:nModel);
Minverted(1:nModel) = [];

MIevt                = Minverted(1:nEvents);
Minverted(1:nEvents) = [];

MIsta   = Minverted;
MImodel = reshape(MImodel,size(xMat));

%% make some pretty plots

figure(1)

if collapse_y

    plot(xMat(1, :), MImodel(1, :));
    xlabel('Distance, km');
    ylabel('\Deltat*')
    
else
    
    contourf(xMat,yMat,MImodel, v); h = colorbar; %caxis([0 2])
    h.Label.String = '\Deltat*';
    load('CMfine');
    colormap(cm)
    caxis(clim)
    hold on

    [dX, indY] = unique(dataX);

    scatter(dX, dataY(indY),10, MIsta,'filled', 'MarkerEdgeColor', 'k')
    set(gca,'YDir','normal')
    daspect([1 1 1])

    C=load('coast');
    [C.x,C.y]=mfwdtran(mstruct,C.lat,C.long);
    C.x(C.x>maxX)=nan;
    C.x(C.x<minX)=nan;
    C.y(C.y>maxY)=nan;
    C.y(C.y<minY)=nan;
    plot(C.x,C.y,'k-','LineWidth',2)

    xlim([ minX maxX ]);
    ylim([ minY maxY ]);

end
    
figure(2)
histogram(MIsta)
xlabel('Size of the station terms');

save(name, 'MImodel', 'xMat', 'yMat', 'mstruct', 'MIevt', 'MIsta', 'dataE', 'dataI', 'dataV', 'dataX', 'dataY', 'uSta');
