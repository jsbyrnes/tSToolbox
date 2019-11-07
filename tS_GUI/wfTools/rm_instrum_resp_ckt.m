%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% rm_instrum_resp_ckt: a script to remove the instrument response for AVO 
% station CKT, a Marks Products L4 seismometer. Telemetry is also accounted
% for in the correction.
%
% Last update 4/22/2012
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% load the raw data
rawdata = load('ckt_20051105022900_20051105023400.txt');

%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% AVO specific information
%
%%%%%%%%%%%%%%%%%%%%%%%%%%

% Damping constant of instrument (underdamped < 1 < overdamped)
beta = 0.67;

% Natural frequency of instrument (Hz)
fnat = 1.0;

% Calibration value (nm/s/count)
calib = 12.303;

% Cutoff frequency of mcvco (Hz)
mcvcof = 30;

% Order of the mcvco
mcvcoord = 2;

% Cutoff frequency of the discriminator (Hz)
discf = 30;

% Order of the discriminator
discord = 2;

% frequency-amplitude pairs for the antialias filter

% amplitudes of generic AVO antialias filter - used at CKT on 11/05/2005
% Type II anti-alias filter
aaresp = [ 1.000 1.000 1.000 0.980 0.980 0.960 0.960 0.960 0.960 0.960 ...
           0.900 0.810 0.720 0.643 0.600 0.551 0.520 0.488 0.456 0.424 ...
           0.392 0.360 0.344 0.328 0.312 0.296 0.280 0.236 0.192 0.148 ...
           0.104 0.060 ];

% amplitudes of J-120D antialias filter - not used at CKT on 11/05/2005
% Type II anti-alias filter
% aaresp = [ 1.000 1.000 1.000 1.000 1.000 1.000 1.000 1.000 1.000 1.000 ...
%            1.000 1.000 1.000 1.000 1.000 1.000 1.000 0.812 0.676 0.549 ... 
%            0.436 0.354 0.288 0.234 0.190 0.154 0.100 0.081 0.068 0.055 ... 
%            0.044 0.029 ];
       
% frequencies of the antialias filter (Hz)
rmat = [   0.100 0.200 0.350 0.500 0.750 1.000 1.670 2.000 2.500 3.330 ...
           5.000 7.500 10.00 12.50 15.00 17.50 20.00 22.00 24.00 26.00 ...
           28.00 30.00 32.00 34.00 36.00 38.00 40.00 42.00 44.00 46.00 ...
           48.00 50.00 ];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Compute analog poles/zeros
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% initialize the size of polz
polz = zeros(2+mcvcoord+discord,1);

% poles for the instrument
% from equation 4.54 on page 64 in the book "Of Poles and Zeros"
% this assumes beta is less than 1 (the underdamped case)
% note that what is called beta here is called h in that book
polz = [ -(beta*2*pi*fnat) + i*sqrt(1-(beta^2))*(2*pi*fnat) ...
         -(beta*2*pi*fnat) - i*sqrt(1-(beta^2))*(2*pi*fnat) ]; 

% poles for the McVCO, from formula for analog Butterworth poles     
for ii=1:mcvcoord
    polz(2+ii) = (mcvcof*2*pi)*exp(i*((2*ii+mcvcoord-1)/(2*mcvcoord))*pi);
end

% poles for the discriminator, from formula for analog Butterworth poles
for ii=1:discord
    polz(2+mcvcoord+ii) = ...
        (discf*2*pi)*exp(i*((2*ii+discord-1)/(2*discord))*pi);
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% 12 pieces of input needed
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%  1. Sample rate (Hz)
samplrate = 100;

%  2. Low frequency cutoff (Hz)
flo = 0.1;

%  3. High frequency cutoff (Hz)
fhi = 10;

%  4. Butterworth filter order at low frequency cutoff (between 2 and 4)
ordl = 3;

%  5. Butterworth filter order at high frequency cutoff (between 3 and 7)
ordh = 5;

%  6. Bad data value (e.g., for telemetry dropouts) - this currently
%  doesn't work for NaN, so you have to change any NaNs to this prior to
%  applying these codes
badvals = (-2^-31);

%  7. Zeros for the instrument (rad/s, not Hz)
zers = [  0 ...
          0 ];
      
%  8. Poles for the instrument (rad/s, not Hz)
pols = polz;
     
%  9. Inverse gain factor (counts/m/s)
digout = calib*(10^-9); 

% 10. Frequency at which normalization of the response occurs (Hz)
digoutf = 5;

% 11. Over-sampling rate for accurate instrument response
ovrsampl = 5;

% 12. Intrinsic delay of the system (ms) - zero for minimum phase, 35 ms 
% for AVO short period stations w/o antialias filter, 55 ms for stations 
% with anialias filter
idelay = 55;


% Instrument correct
prcdata = rm_instrum_resp(rawdata,badvals,samplrate,pols,zers,...
    flo,fhi,ordl,ordh,digout,digoutf,ovrsampl,idelay);


% Account for the antialias filter: create an approximate digital filter  
hl = firls(13-1,rmat/(samplrate/2),aaresp);
[a, b] = freqz(hl,1,samplrate);

% Compute the minimum phase equivalent of the digital filter
mp = real(ifft(abs(fft(hl)).*exp(-i*imag(hilbert(log(abs(fft(hl))))))));

% Apply the minimum phase filter
dum = conv(mp,prcdata);
prcdata = dum(1:length(prcdata));


% Plot the instrument corrected trace
plot([0:(length(prcdata(1:end))-1)]*(1/samplrate),...
    prcdata(1:end)*(10^6),'r','LineWidth',2)











