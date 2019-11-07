function [E] = findEvents(srchParams)

%

%unpack search parameters and use defaults as needed.
sTime=srchParams.sTime;
eTime=srchParams.eTime;
radCoords=srchParams.radialCoords;
minZ=srchParams.minZ;
maxZ=srchParams.maxZ;
minM=srchParams.minM;
maxM=srchParams.maxM;

E=irisFetch.Events('radialcoordinates'...
    ,radCoords,'startTime',sTime,'endTime',eTime,'minimumMagnitude',minM,...
    'maximumMagnitude',maxM,'minimumDepth',minZ,'maximumDepth',maxZ);



