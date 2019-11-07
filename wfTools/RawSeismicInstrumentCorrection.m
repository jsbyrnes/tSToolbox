function [ correctedTrace polesVector zerosVector ] = RawSeismicInstrumentCorrection(inputTrace, polesVector, zerosVector, polesZerosUnits, sampleRate, minFreq, maxFreq, sensitivity, plotOverview)
% RawSeismicInstrumentCorrection
% (inputTrace,polesVector, zerosVector, polesZerosUnits, sampleRate, minFreq, maxFreq ,sensitivity, plotOverview)
% corrects raw record for amplitude and phase characteristics of the instrument.
%
% Correction is applied in spectral domain
% to frequency band from desired minFreq
% and MaxFreq or Nyquist frequency if MaxFreq > Nyquist frequency.
%
% input : input trace
%         vector of poles
%         vector of zeros 
%         poles and zeros units - string - must be set to either 'Hz' or 'rad'
%         sampleFrequency [Hz]
%         minimal frequency threshold for transfer function to apply
%         maximal frequency threshold for transfer function to apply
%         sensitivity constant = sensor constant * digitizer constant
%         plotOverview - type 'plot' to see overview figure or something
%                        else to disable results plot
%         
% output: vector of corrected signal
%         vector of poles [Hz]
%         vector of zeros [Hz]
%
% Important Note:
% To work properly the prior filtration of the input trace for the desired 
% frequency range is required. I have decided to let the user select 
% the filter which fits best her/his needs. 
%
% 
%
% credit: Martin Mityska (2014)

    signalLength=numel(inputTrace);

    % Let's remove mean
    inputTrace=inputTrace-mean(inputTrace);

    % Let's remove linear trend
    inputTrace=detrend(inputTrace);

    % Step in frequency domain
    deltaFreq = sampleRate/(signalLength);
    
    nyquistFreq = sampleRate/2; 

    % Is maxFreq lower than Nyquist frequency?
    if maxFreq > nyquistFreq
        warning(...
        ['Correction maxFreq is larger than Nyquist frequency of the signal.\n'...
         'Requested maxFreq: %f, Nyquist freq: %f\n'...
         'maxFreq set to Nyquist frequency.']...
          ,maxFreq,nyquistFreq)
      maxFreq = nyquistFreq;
    end
    
    % Obtaining transfer function in frequency domain for given zeros, poles, dF and frequency band
    [transfer polesVector zerosVector]=GetInstrumentTransferFunc(polesVector, zerosVector,polesZerosUnits, minFreq ,maxFreq, deltaFreq,sensitivity);

    transfer=[ones(round(minFreq/(deltaFreq)-1),1).*1.0/sensitivity;transfer];

    % Let's get spectrum of the signal given
    spectrum=fft(double(inputTrace));
    
    transferFitted=[];
    if numel(transfer) ~= numel(spectrum)/2
        misshift = numel(transfer)-numel(spectrum)/2;
        if misshift < 0
            transferFitted=[transfer;ones(abs(misshift),1).*1.0/sensitivity];
        end
        if misshift > 0
            transferFitted=[transfer(1:numel(transfer)-misshift)];

        end
    else
       transferFitted=transfer;
    end

    % Building shape of transfer function - adding redundant part needed 
    % for dividing of the spectrum and setting imag part to odd function
    % (necessary for correct phase correction).
    %
    % Test for odd length of the input trace.
    if(mod(signalLength,2)~=0)
        % ->                                                                   %-Fix for odd length of input trace ---------------------------%             
        transferFitted=[complex(real(transferFitted),-1*imag(transferFitted)); complex(real(transferFitted(end)),-1*imag(transferFitted(end))); flipud(complex(real(transferFitted),imag(transferFitted)))];
    else
        % If the length of the trace is even, the spectrum and transfer
        % function is symmetric and no fix is needed.
        transferFitted=[complex(real(transferFitted),-1*imag(transferFitted));flipud(complex(real(transferFitted),imag(transferFitted)))];
    end
    
    % In this step, we divide spectrum of the signal by the transfer
    % function. The integration in frequency domain could be done now. It
    % depends actually on the given number of zeros.
    spectrumCorrected = spectrum./transferFitted;
    
    % The signal corrected for instrument response in time domain is
    % obtained by application of inversion Fourier transform. The real part
    % of the result represents the requested signal.
    correctedTrace = real(ifft(spectrumCorrected));
    
    if(strcmpi(plotOverview,'plot'))
        PlotInstrumentTransferFunc(polesVector, zerosVector, polesZerosUnits, minFreq, maxFreq, deltaFreq, sampleRate, sensitivity, inputTrace, correctedTrace,transferFitted)
    end
    
end

function [ transfer, polesVector, zerosVector ] = GetInstrumentTransferFunc(polesVector, zerosVector, polesZerosUnits, minFreq ,maxFreq, deltaFreq, sensitivity)
% GetInstrumentTransferFunc
%(polesVector, zerosVector, polesZerosUnits, minFreq ,maxFreq, deltaFreq,sensitivity)
%
% returns complex transfer function for selected
% instrument and for frequency band from minF to
% maxF with frequency step deltaF.
%
%         vector of poles
%         vector of zeros 
%         poles and zeros units - string - must be set to either 'Hz' or 'rad'
%         minFrequency             [Hz]
%         maxFrequency             [Hz]
%         step in frequency domain [Hz] 
%
% output: complex transfer function
%         instrument database entry - struct
%         vector of poles          [Hz]
%         vector of zeros          [Hz]
%
% The transfer function in frequency domain is calculated for given poles,
% zeros and sensitivity constant.
%                                 ---
%                                 | |Nz   
%                 1          -1 * | |z=1 ( if - Zz )
% H(f[Hz]) = -----------  * --------------------------
%            sensitivity          --- 
%                                 | |Np   
%                                 | |p=1 ( if - Pp )
% 
% So, the transfer function in frequency domain is calculated as -1* product of
% roots of zeros polynomial divided by product of roots of poles polynomial.
% Poles and zeros are roots of the polynomials.
% 
% Reference:
% Frank Scherbaum. 2001: Of Poles and Zeros: Fundamentals of Digital Seismology (2nd ed.).
% Kluwer Academic Publishers, Norwell, MA, USA.
%
% credit: Martin Mityska (2014)


	% We have to check units of poles and zeros
    if(strcmpi(polesZerosUnits,'rad') )
        
        % Let's convert poles and zeros from radians to Hz
        polesVector = polesVector./(2*pi());
        zerosVector = zerosVector./(2*pi());
        
        else if((strcmpi(polesZerosUnits,'Hz') ))
            % Poles and zeros are already in Hz, nothing has to be done.
                1;
        else
            
        error('Undefined poles and zeros dimension - the 3rd argument has to be set to either ''rad'' or ''Hz'' ');
    
        end
    end
    
    fAxis = [minFreq:deltaFreq:maxFreq];

    % Here we construct numerator and denominator of the transfer equation:
        numerator=-1.*ones(1,numel(fAxis));
        denominator=ones(1,numel(fAxis));

        for k=1:numel(zerosVector)
            numerator=numerator.*( fAxis.*1i - zerosVector(k) );
        end

        for k=1:numel(polesVector)
            denominator = denominator.*( fAxis.*1i - polesVector(k) );
        end

    transfer = ((numerator./denominator)')/sensitivity;

end

function [hFig ] = PlotInstrumentTransferFunc(polesVector, zerosVector, polesZerosUnits, minFreq, maxFreq, deltaFreq, sampleRate, sensitivity, originalTrace, correctedTrace, transferUsed )
% PlotInstrumentTransferFunc
% (polesVector, zerosVector, polesZerosUnits, minFreq, maxFreq, deltaFreq, sampleRate, sensitivity, originalTrace, correctedTrace, transferUsed )
% plots results overview for selected instrument correction
% to loglog figure.
%
% input : station code - string
%         channel code - string
%         minFrequency             [Hz]
%         maxFrequency             [Hz]
%         step in frequency domain [Hz] 
% output: figure handle
%
% credit: Martin Mityska (2014)

    hFig=figure();

    set(hFig, 'Position', [0 0 1024 768]);
    set(gcf,'color','w');

    minFreqForInformativeDisplay = 0.01; %Hz
    
    % Transfer function is displayed for overview only - range
    % minFreqForInformativeDisplay - Nyquist frequency
    [transfer polesVector zerosVector]=GetInstrumentTransferFunc(polesVector, zerosVector,polesZerosUnits, minFreqForInformativeDisplay ,sampleRate/2, deltaFreq,sensitivity);

    dt = 1.0/deltaFreq;
    
    fAxis = [minFreqForInformativeDisplay:deltaFreq:sampleRate/2];
    tAxis = [(0:1:numel(originalTrace)-1).*dt];
    
    subplot(5,3,1);
   
    loglog((fAxis),(abs(transfer)),'LineWidth',2);
    title(sprintf('Instrument transfer function\n'),'fontsize', 16);
    xlabel([' f [Hz] ']);
    ylabel('abs(transfer function)');
    set(gca,'XTick',[10^-3,10^-2,10^-1,10^0,10^1,10^2,10^3]);
    
    set(gca,'XMinorTick','on','YMinorTick','on')
    grid on;

    subplot(5,3,4);
    
    phase = unwrap(atan2(imag(transfer),real(transfer)));
    semilogx(fAxis,phase,'LineWidth',2);
    title(sprintf('Phase characteristic of the instrument\n'),'fontsize', 16);
    xlabel([' f [Hz]  ']);
    ylabel('phase [rad]');
    set(gca,'XTick',[10^-3,10^-2,10^-1,10^0,10^1,10^2,10^3]);
    
    set(gca,'XMinorTick','on','YMinorTick','on')
    grid on;
    
    subplot(5,3,2:3)
    
    plot(tAxis,originalTrace);
    title('Input trace','fontsize', 16);
  
    
    subplot(5,3,5:6)
    hold on
    plot(tAxis,correctedTrace,'b');
    hold off
    title('Corrected trace','fontsize', 16);
    
    subplot(5,3,8:9)
    hold on
    plot(tAxis,normalize(originalTrace,1),'b');
    plot(tAxis,normalize(correctedTrace,1),'r');
    hold off
    title('Corrected (red) and original trace (blue) - normalized','fontsize', 16);
    
    subplot(5,3,11:12)
    
    semilogy((0:1:numel(transferUsed)-1)*deltaFreq,abs(transferUsed))
    xlabel([' f [Hz]  ']);
    title('Absolute value of actually used transfer function (+ redundant part)')
    
    subplot(5,3,14:15)
    
    plot((0:1:numel(transferUsed)-1)*deltaFreq,imag(transferUsed))
    xlabel([' f [Hz]  ']);
    title('Imaginary part of actually used transfer function (+ redundant part)')
    set(gcf,'Color',[1 1 1]);
    
    h= legend(PrepareLegendText(polesVector, zerosVector, minFreq, maxFreq,sensitivity, sampleRate));
 
    legendPosRect = [0.1, 0.25, .25, .1];
    set(h, 'Position', legendPosRect)
    set(h,'fontname', 'Courier');
    set(h,'fontweight', 'Bold');
    set(h,'fontsize', 10);

end

function [legText] = PrepareLegendText(pole, zero, minFreq, maxFreq, sensitivity, sampleRate)
    
    legText = [sprintf('sample frequency: \n')];
    legText = [legText sprintf('%10.5f Hz\n', sampleRate)];
    legText = [legText sprintf('Frequency range\nfor correction:\n')];
    legText = [legText sprintf('%9.5f  -  %9.5f Hz\n',minFreq,maxFreq)];

        legText = sprintf([legText '\n poles     :\n']);
    
    for k = 1:numel(pole)
        if real(pole(k)) < 0.0;
            legText = [legText sprintf('%010.5f ',real(pole(k)))];
        else
            legText = [legText sprintf('%09.5f ',real(pole(k)))];
        end
        if imag(pole(k)) < 0.0
            legText = [legText sprintf('- %09.5fi\n',abs(imag(pole(k))))];
        else
            legText = [legText sprintf('+ %09.5fi\n',imag(pole(k)))];
        end
    end
    
    legText=[legText sprintf(' zeros     :\n')];
    
    
    for k = 1:numel(zero)
        if real(zero(k)) < 0.0
            legText = [legText sprintf('%10.5f ',real(zero(k)))];
        else
            legText = [legText sprintf(' %09.5f ',real(zero(k)))];
        end
        if imag(zero(k)) < 0.0
            legText = [legText sprintf('- %09.5fi\n',abs(imag(zero(k))))];
        else
            legText = [legText sprintf('+ %09.5fi\n',imag(zero(k)))];
        end
    end
    
    legText = [legText sprintf('\nsensitivity: \n')];
    legText = [legText sprintf('%10.5E \n', sensitivity)];  
    
end

function [normTracks] =  normalize(data,params)
% Normalization function
%
% input:
%        data   - Column major matrix of traces.
%        params - 0 = Normalization according to global max value of all traces.
%               - 1 = Normalization of every single trace separately.
%
% Martin Mityska (2011)

    normTracks=zeros(size(data,1),size(data,2));

    % Normalization of all traces according to max value of all of them.
    if params == 0;
          maxValueOfTracks = max(max(abs(data)));
          for trackIndex=1:1:length(data(1,:));
              if maxValueOfTracks ~= 0.0 ;
                 normTracks(:,trackIndex)=(data(:,trackIndex)/maxValueOfTracks);
              end
          end
    end

    % Applying normalization to every single trace separately.
    if params == 1 
          maxValueInTrack = max(abs(data));
          for trackIndex=1:1:length(data(1,:));
                normTracks(:,trackIndex)=(data(:,trackIndex)/maxValueInTrack(trackIndex));
          end
    end
    
end		