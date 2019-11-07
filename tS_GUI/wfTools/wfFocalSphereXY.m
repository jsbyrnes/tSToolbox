function [ wf ] = wfFocalSphereXY( wf, evt )
% wfFocalSphereXY, for a given set of (irisFetch) waveforms and an event
% structure, this function calculates the X and Y position on the focal
% sphere (lower hemisphere) of each recording.

load myTakeoffAngles

d=round(evt.depth);

%here take just a column of the takeoff matrix, corresponding to the depth
%of the event (rounded to the nearest km).
takeOff=tkof(:,depth==d);
%this is the delta corresponding to the transition from upgoing to
%downgoing  
cutoff=cutoff(depth==d);

for k=1:length(wf)
    
    %Del is the epicentral distance for this station
   [Del azm]=distance(evt.lat,evt.lon,wf(k).latitude,wf(k).longitude);
   depth=evt.depth;
   
   tkAngle=interp1(delta,takeOff,Del); 
   
    x = sind(azm)*sind(tkAngle); 
    y = cosd(azm)*sind(tkAngle); 
    z = cosd(tkAngle); z=-z;
    
    % if rays are going up, I have to filp the sign of x and y because I'm
    % porjecting onto lower hemisphere
    if Del < cutoff; 
        x=-x; y=-y; 
    end
    
    
    px=x*sqrt( 2/(1-z) ); py=y*sqrt( 2/(1-z) );
   
    wf(k).focalSphereXY=[px py];
    
    
    
end




