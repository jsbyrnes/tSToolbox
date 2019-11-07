%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% rm_instrum_resp_spcr: a script to remove the instrument response for AVO 
% station SPCR, a Guralp 6TD seismometer.
%
% Last update 4/22/2012
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% load the raw data
rawdata = load('spcr_20051105022900_20051105023400.txt');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% 12 pieces of input needed
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%  1. Sample rate (Hz)
samplrate = 50;

%  2. Low frequency cutoff (Hz)
flo = 0.1;

%  3. High frequency cutoff (Hz)
fhi = 10;

%  4. Butterworth filter order at low frequency cutoff (between 2 and 4)
ordl = 3;

%  5. Butterworth filter order at high frequency cutoff (between 3 and 7)
ordh = 5;

%  6. Bad data value (e.g., for telemetry dropouts) - this doesn't 
%  work for NaN, so you have to change any NaNs to this prior to
%  applying these codes
badvals = (-2^-31);

%  7. Zeros for the instrument (rad/s, not Hz)
zers = 2*pi*...
       [ -5.03207 ...
          0 ...
          0 ];
      
%  8. Poles for the instrument (rad/s, not Hz)
pols = 2*pi*...
       [ -23.65*(10^-3) + i*23.65*(10^-3) ...
         -23.65*(10^-3) - i*23.65*(10^-3) ...
         -393.011 ...
         -7.4904 ...
         -53.5979 - i*21.7494 ...
         -53.5979 + i*21.7494 ];
     
%  9. Inverse gain factor (m/s/count)
digout = 2.326*(10^-10);

% 10. Frequency at which normalization of the response occurs (Hz)
digoutf = 1;

% 11. Over-sampling rate for accurate instrument response
ovrsampl = 5;

% 12. Intrinsic delay of the system (ms) - zero for minimum phase, 35 ms 
% for AVO short period stations w/o antialias filter, 55 ms for stations 
% with anialias filter 
idelay = 0;

% Instrument correct
prcdata = rm_instrum_resp(rawdata,badvals,samplrate,pols,zers,...
    flo,fhi,ordl,ordh,digout,digoutf,ovrsampl,idelay);

% make a figure with maximum screen extent
scrsz = get(0,'ScreenSize');
figure('Position',[1 scrsz(1) scrsz(3) scrsz(4)])

fsize = 24;
plot([0:(length(prcdata)-1)]*(1/samplrate),prcdata*(10^6),'LineWidth',2)
axis([0 300 -6*(10^0) 6*(10^0)]); hold on
set(gca,'Fontsize',fsize,'FontWeight','bold');
xlabel(' Time (s) '); ylabel(' Particle velocity (\mum/s) '); 
title(' Broadband (blue)/Short-period (red) Comparison ')

