function  [TrCorr] = wfInstCorr1trace_acausal(Tr, resp)

   %zero-pole to transfer function conversion (MATLAB built-in)
   [resp.poly_num, resp.poly_den] = zp2tf(resp.Zeros*2*pi, resp.Poles*2*pi, resp.Amp); 
   % matlab function ADDED 2PI because zeroes and poles are in Hz in RESP files, need them in rad/s
   
   Z=Tr.data;
   
   npts = length(Z);
   bfwave=zeros(1,npts);
   
   Zfftdata=fft(Z);
   fftlength=npts;
   
   dt = 1/Tr.sampleRate;
   df = 1/dt;

   % remove instrument response for the waves spectrum=zeros(1,npts);

    f = df*(1:fftlength)/fftlength;  % half length
    w = 2*pi*f'; % angular freq. w=2*pi*f; I need for this to be a column
        
    h = freqs(resp.poly_num, resp.poly_den, w); % frequency response (built-in)

    % remove instrument response: the first half frequency spectrum
    Zfftdata = Zfftdata.*conj(h)./max(abs(h).^2, 0.01);  
       
    Zcorr = real(ifft(Zfftdata));

    TrCorr=Tr;
    TrCorr.data=Zcorr;
