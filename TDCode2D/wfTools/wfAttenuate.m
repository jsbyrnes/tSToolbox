function [ Ta ] = wfAttenuate( T, tStar, w1)
% this function takes an irisFetch trace structure (or array of structures)
% and returns a similar (array of)
% structure where the 'data' has been attenuated by a t* indicated as input
%
% USAGE 
% [ Ta ] = wfFFT( T, tStar )
% where T is a trace structure.
% Ta for T-attenuated.
% it will add a field to the structure called tStar to hold the value used

%first check to see if a spectrum has already been calculated,
%otherwise, calculate it

if ~isfield(T,'ampSpec')
    T=wfFFT2(T);
end

Ta=T;
for k=1:length(T)
        
    %Calculate the freq-domain attenuation operator:
    
    if nargin ~= 3
    
        w1 = 1*2*pi; % reference frequency (at unrelaxed modulus), 400 set before. 0.1 gives smoothest curve
    
    end
    
    % set DC term to small non-zero number (divide protect)
    T(k).ampSpectrum(1)=T(k).ampSpectrum(2)*.001;
    
    o = 0; %offset that goes into attn, I'm leaving as 0
    
    attnOperator = calcAttnOperator(tStar,2*pi*T(k).frequencies,w1,o);
    
    Ta(k).data=real( ifft( fft(T(k).data).*attnOperator ) );
    Ta(k).tStar=tStar;
    
    %recalculate Fourier transform
    Ta(k)=wfFFT2(Ta(k));
        
end

function attnOperator = calcAttnOperator(tStar,frequencies,w1,offset)
%
% [AttF] = attn(t_star,W,w1,offset)
%  AttF is freq-domain attenuation opperator
%  t_star is t*
%  W is array of circular frequencies where calculations are made
%    (from zero to Nyq freq, higher freqs are made by reflection)
%  w1 is reference frequency
%  offset is amount of offset of the atttenuation pulse

%JSB alpha is not yet in use

 if nargin==2       % set w1 AND offset to default values
     w1=400;
     offset=0;
 end
 
 if nargin==3       % set w1 OR  offset to default values
   if w1>10 % case where w1 is w1
     offset=0;
   else % case where w1 is offset
     offset=w1;
     w1=400*2*pi;
   end
 end


%-------------------------------------------------------------------
Amp = exp(-frequencies.*tStar/2);                 % amplitude spectrum
Dt  = -(tStar/pi).*log([1,frequencies(2:end)]/w1); % dummy F=0, OK on next line
Dph = frequencies.*(Dt)+offset;                   % delta phase shift spectrum

attnOperator= Amp.*(cos(Dph) + 1i*sin(Dph));     % attn opperator in freq
attnOperator=attnOperator';
%-------------------------------------------------------------------


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%% Notes %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% t*star is the input parameter
%
%%%%%% deriving Amp and (in two ways) dispersion relation,
%    (with relations that assume freq independent Q)
%
% 1) using c/c1 = 1 + (1/(pi*Q))*ln(W/w1)
%    t_star=x/(c*Q) => Q=x/(c*t_star)
%    alpha = W/(2*c*Q) = W*t_star/(2*x);
%    Amplitude is exp(-alpha*x) = exp(-W*t_star/2)
%    c = c1*{1+ln(W/w1)/(pi*Q)}
%    t = x/c = x/{c1*[1+ln(W/w1)/(pi*Q)]}
%    t  .appx. (x/c1)*{1-ln(W/w1)/(pi*Q)}
%    dt .appx. -[x/(c1*pi*Q)]*ln(W/w1)
%    dt .appx. -(t_star/pi)*ln(W/w1)     <========
%    dPhi = -dt*W
%
% 2) Azumi has alpha = a0*W  and H{alpha} = -[2*a0*W/pi]*ln(a1*W)
%    since alpha=W/(2*c*Q), a0=1/(2*c*Q) = t_star/(2*x)
%    Then, H{alpha} = -t_star*W/(pi*x)*ln(a1*W)
%    (1/a1 is recognized as the  hi-cut corner)
%    Using K-K, c = 1/(1/c_inf + H{alpha}/W)
%    t = x/c = x*{1/c_inf - [t_star/(pi*x)]*ln(a1*W}
%    t = x/c_inf - (t_star/pi)*ln(a1*W)
%    dt = -(t_star/pi)*ln(a1*W)           <========
%    dPhi = -dt*W

