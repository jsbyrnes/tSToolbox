function  [TrCorr] = wfInstCorr1trace_causal(Tr, resp)

   %zero-pole to transfer function conversion (MATLAB built-in)
   [resp.poly_num, resp.poly_den] = zp2tf(resp.Zeros*2*pi, resp.Poles*2*pi, resp.Amp); 
   % matlab function ADDED 2PI because zeroes and poles are in Hz in RESP files, need them in rad/s
   
   TrCorr = Tr;
   
    TrCorr.data = rm_instrum_resp(Tr.data,1e10,Tr.sampleRate,resp.Poles*2*pi,resp.Zeros*2*pi,...
        1/200,10,3,4,1/Tr.sensitivity,Tr.sensitivityFrequency,1,0);
