%% NMR_TimeFit
% A class that fits FID to a series of exponentially decaying components using a
% NMR_Mix object. This class knows how to fit to a time domain signal with
% or without constraints.
%
classdef NMR_TimeFit < NMR_Fit
    methods
        function obj = NMR_TimeFit(time_domain_signal, t, ...
                area, freq, fwhm, phase, line_broadening, zeroPadSize)
            % Construct an NMR_TimeFit
            obj = obj@NMR_Fit(time_domain_signal, t, ...
                area, freq, fwhm, phase, line_broadening, zeroPadSize);
        end
        function [fit_area, fit_freq, fit_fwhm, fit_phase] = ...
                fitComponentToResidual(obj)
            % Fits a single exponentially decaying component to the
            % residual signal. This is useful in peak fitting. Note,
            % however, that this function only fits based on residuals,
            % thus you are fitting the residuals, and not the signal. After
            % fitting the residuals, it is recommended that you "refit" all
            % the components as they will likely change slightly to better
            % accomodate the new component.
            
            %% Step one:
            % Use the magnitude of the residual freq domain signal to
            % make a guess of the next peaks frequency
            
            % If there are already components, fit them first, then
            % subtract them to get the residual spectrum
            if(~isempty(obj.area))
                residualSpectrum = obj.spectralDomainSignal - ...
                    obj.calcSpectralDomainSignal(obj.f);
            else
                residualSpectrum = obj.spectralDomainSignal;
            end
            
            % Find a guess frequency for the next peak
            [maxVal maxIdx] = max(abs(residualSpectrum));
            peakFreq = obj.f(maxIdx);
            
            %% Step two:
            % Fit the residual signal in the time domain, giving reasonable
            % bounds
            freq_halfBound = 100;
            fwhm_UpperBound = 500;
            residualTimeDomainSignal = obj.timeDomainSignal - obj.calcTimeDomainSignal(obj.t);
            largestPeakFit = NMR_TimeFit(residualTimeDomainSignal, obj.t, ...
                1, peakFreq, 100, 0, obj.lineBroadening, obj.zeroPadSize);
            largestPeakFit = largestPeakFit.setBounds(0, inf, ...
                peakFreq-freq_halfBound, peakFreq+freq_halfBound,...
                0, fwhm_UpperBound,-inf,inf);
            [fit_area, fit_freq, fit_fwhm, fit_phase] = ...
                largestPeakFit.calcTimeDomainSignalFit();
        end
        
        function obj = autoAddComponent(obj)
            % This function attemps to automatically fit and add a new
            % exponentially decaying signal to the NMR_mix object
            
            % Fit next component to residual signal, then add it to this
            % NMR_Mix object
            [add_area, add_freq, add_fwhm, add_phase] = ...
                obj.fitComponentToResidual();
            
            % Add fitted component to NMR_mix
            
            obj = obj.addComponents(add_area, add_freq, add_fwhm, add_phase);
            
            % Refit all components after addition of latest component
            obj = obj.fitTimeDomainSignal();
            
            if(~isempty(obj.ub))
                % Remove bounds in case sorted order changed
                obj.ub = [];
                obj.lb = [];
                warning('Deleting bounds - you need to call setBounds again with the new components bounds');
            end
        end
        
        function obj = autoAddComponents(obj, nComponents)
            % This function attemps to automatically fit and add
            % n (defined by nComponents) new exponentially decaying signals
            % to the NMR_mix object
            
            for iComp = 1:nComponents
                obj = obj.autoAddComponent();
            end
        end
        
        function obj = fitTool(obj)
            % This function attemps to adds exponentially decaying signal
            % components until the user stops the fit.
            
            % Create a figure to show results
            thisFigure = gcf();
            clf;
            
            % Keep fitting until user says stop
            continueFitting = true;
            while(continueFitting)
                % Create temporary NRM_Fix object and add newest component
                tempFit = obj;
                tempFit = tempFit.autoAddComponent();
                
                % Show fit
                tempFit.plotFit();
                
                % Report new fits
                tempFit.describe();
                
                % If user wants to keep component, add it to the object
                output_str1 = lower(input('Keep fit? [y/n]','s'));
                keepFit = length(findstr(output_str1, 'y')>0);
                if(keepFit)
                    obj = tempFit;
                end
                
                % If user wants to fit more peaks, continue fitting...
                output_str2 = lower(input('\tFit more peaks? [y/n]','s'));
                continueFitting = length(findstr(output_str2, 'y')>0);
            end
        end
        
        function [fit_area, fit_freq, fit_fwhm, fit_phase] = calcTimeDomainSignalFit(obj)
            % Fits exponentially decaying components to the given time
            % domain signal and provides the results without saving them to
            % this object. To save the results to the NMR_Mix, use
            % fitTimeDomainSignal
            
            fitoptions = optimoptions('lsqcurvefit');
            %                         fitoptions.Display = 'iter-detailed';
            %                         fitoptions.Display = 'final-detailed';
            fitoptions.Display = 'off';
            fitoptions.MaxIter = 10000;
            fitoptions.TolFun=1E-900;
            fitoptions.TolX = 1E-15;
            fitoptions.FinDiffType = 'central';
            fitoptions.Algorithm = 'trust-region-reflective';
            fitoptions.MaxFunEvals = 5000;
            if(isempty(obj.lb) & isempty(obj.ub))
                % We can use complex fitting if we dont want constraints...
                % this should go faster
                
                % Put all components into a single matrix
                guess = [obj.area.*exp(1i*pi*obj.phase/180); 1i*2*pi*obj.freq-pi*obj.fwhm];
                                
                fit_params = lsqcurvefit(@obj.calcUnconstrainedTimeSig,guess,obj.t,...
                    obj.timeDomainSignal,[],[],fitoptions);
                
                % Separate out the components from the matrix
                fit_area = abs(fit_params(1,:));
                fit_freq = imag(fit_params(2,:))/(2*pi);
                fit_fwhm = -real(fit_params(2,:))/pi;
                fit_phase = angle(fit_params(1,:))*180/pi;
            else
                % Put all components into a single matrix
                guess = [obj.area; obj.freq; ...
                    obj.fwhm; obj.phase];
                
                fit_params = lsqcurvefit(@obj.calcConstrainedTimeSig,guess,obj.t,...
                    [real(obj.timeDomainSignal),imag(obj.timeDomainSignal)],...
                    obj.lb,obj.ub,fitoptions);
                
                % Separate out the components from the matrix
                fit_vec = fit_params(1,:).*exp(1i*pi*fit_params(4,:)/180);
                fit_area = abs(fit_vec);
                fit_freq = fit_params(2,:);
                fit_fwhm = fit_params(3,:);
                fit_phase = angle(fit_vec)*180/pi;
            end
            
            
        end
        
        function obj = fitTimeDomainSignal(obj)
            % Fits exponentially decaying components to the given time
            % domain signal and saves the results to this NMR_Mix object.
            % To just return the fits and not save the results, use
            % calcTimeDomainSignalFit
            [fit_area, fit_freq, fit_fwhm, fit_phase] = obj.calcTimeDomainSignalFit();
            
            % Save fits
            obj = obj.resetComponents(fit_area, fit_freq, fit_fwhm, fit_phase);
        end
        
        function complexSig = calcUnconstrainedTimeSig(obj,nmr_params,t)
            % A function used in fitting to allow constraints for complex
            % fitting. This is the same as calling calcTimeDomainSignal of
            % the NMR_Mix object, only it separates the real and imaginary
            % parts to allow for constraints to be used in fitting.
            nComp = size(nmr_params,2);
            complexSig = zeros(size(t));
            for iComp = 1:nComp
                complexSig = complexSig + nmr_params(1,iComp).*exp(nmr_params(2,iComp)*t);
            end
        end
        
        function realImagSig = calcConstrainedTimeSig(obj,nmr_params,t)
            % A function used in fitting to allow constraints for complex
            % fitting. This is the same as calling calcTimeDomainSignal of
            % the NMR_Mix object, only it separates the real and imaginary
            % parts to allow for constraints to be used in fitting.
            nComp = numel(nmr_params)/4;
            nmr_params = reshape(nmr_params,[4 nComp]);
            tmpNmrMix = NMR_Mix(nmr_params(1,:), nmr_params(2,:), ...
                nmr_params(3,:), nmr_params(4,:));
            complexSig = tmpNmrMix.calcTimeDomainSignal(t);
            realImagSig = [real(complexSig) imag(complexSig)];
        end
        
        function ax1 = plotFit(obj)
            % Calculate fitted and residual spectrums usign time domain
            % signal so that amplitudes are correct even if signal is
            % truncated at the end
            dwell_time = (obj.t(2)-obj.t(1));
            zeroPaddedTime = min(obj.t(:)) + dwell_time*((1:obj.zeroPadSize)-1)';
            
             % Calculate spectrum
            zeroPaddedFreq = linspace(-0.5,0.5,obj.zeroPadSize+1)/dwell_time;
            zeroPaddedFreq = zeroPaddedFreq(1:(end-1)); % Take off last sample to have nSamples  
            zeroPaddedFreq = zeroPaddedFreq(:); % Make a column vector
            
            fittedSpectrum = dwell_time*fftshift(fft(obj.calcTimeDomainSignal(zeroPaddedTime)));
            individualSpectrums = dwell_time*fftshift(fft(obj.calcComponentTimeDomainSignal(zeroPaddedTime),[],1),1);
%             individualSpectrums = obj.calcComponentSpectralDomainSignal(zeroPaddedFreq);
%             fittedSpectrum = obj.calcSpectralDomainSignal(zeroPaddedFreq);
            residualSpectrum = obj.spectralDomainSignal - fittedSpectrum;
            
            % Calculate lorentzian curves for each component
            nComponents = length(obj.area);
            fMat = repmat(zeroPaddedFreq,[1 nComponents]);
            
            legendStrings = cell(1, nComponents);
            for iComp=1:nComponents
                legendStrings{iComp} = ['C\_' sprintf('%03.0f',iComp)];
            end
            
            % Show results to user
            ax2 = subplot(4,1,2);
            plot(zeroPaddedFreq,abs(obj.spectralDomainSignal),'-b');
            hold on;
            plot(zeroPaddedFreq,abs(fittedSpectrum),'-g');
            plot(zeroPaddedFreq,abs(residualSpectrum),'-r');
            hold off;
            ylabel('Magnitude Intensity');
            set(ax2,'xticklabel',{[]}) ;
            set(ax2,'XDir','reverse');
            
%             ax3 = subplot(5,1,3);
%             phaseSig = angle(obj.spectralDomainSignal);
%             phaseFit = angle(fittedSpectrum);
%             phaseDelta = phaseSig - phaseFit;
%             plot(zeroPaddedFreq, unwrap(phaseSig),'-b');
%             hold on;
%             plot(zeroPaddedFreq,unwrap(phaseFit),'-g');
%             plot(zeroPaddedFreq,unwrap(phaseDelta),'-r');
%             plot(zeroPaddedFreq,unwrap(angle(residualSpectrum)),'-c');
%             
%             hold off;
%             ylabel('Phase (Radians)');
%             set(ax3,'xticklabel',{[]}) ;
%             set(ax3,'XDir','reverse');
            
            ax4 = subplot(4,1,3);
            plot(zeroPaddedFreq,real(obj.spectralDomainSignal),'-b');
            hold on;
            plot(zeroPaddedFreq,real(fittedSpectrum),'-g');
            plot(zeroPaddedFreq,real(residualSpectrum),'-r');
            hold off;
            ylabel('Real Intensity');
            set(ax4,'xticklabel',{[]});
            set(ax4,'XDir','reverse');
            
            ax5 = subplot(4,1,4);
            plot(zeroPaddedFreq,imag(obj.spectralDomainSignal),'-b');
            hold on;
            plot(zeroPaddedFreq,imag(fittedSpectrum),'-g');
            plot(zeroPaddedFreq,imag(residualSpectrum),'-r');
            hold off;
            xlabel('Spectral Frequency (Hz)');
            ylabel('Imaginary Intensity');
            legend('Measured','Fitted','Residual');
            set(ax5,'XDir','reverse');
            
            %             if(~isempty(obj.fref))
            %                 % Add PPM axis
            %                 ax1ppm = subplot(5,1,1);
            %                 set(ax1ppm,'units','normalized',...
            %                     'XAxisLocation','top','YAxisLocation','right',...
            %                     'YTick',[],'YTickLabel',[],'Color','none');
            %
            %                 ax1 = axes('Position',get(ax1ppm,'Position'));
            %             else
            ax1 = subplot(4,1,1);
            %             end
            
            plot(fMat,real(individualSpectrums));
            legend(legendStrings);
            ylabel('Component Intensity');
            set(ax1,'xticklabel',{[]}) ;
            set(ax1,'XDir','reverse');
            %             if(~isempty(obj.fref))
            %                 set(ax1ppm,'XDir','reverse');
            %             end
            
            % Keep all x axes in sinc
            linkaxes([ax1,ax2,ax4 ax5],'x');
            
            %             if(~isempty(obj.fref))
            %                 % Initialize XLim in correct units
            %                 set(ax1ppm,'xlim',NMR_Mix.relFreqToPpm(get(ax2,'XLim'),obj.fref));
            %
            %                 % Keep ppm axiz in sinc with freq axis
            %                 xLimListener = addlistener( ax1, 'XLim', 'PostSet', ...
            %                     @(src,evt) set(ax1ppm,'XLim',...
            %                     NMR_Mix.relFreqToPpm(get(ax1,'XLim'),obj.fref)) );
            %             end
        end
        
        function ax1 = plotTimeFit(obj)
            % Calculate fitted and residual spectrums
            dwell_time = (obj.t(2)-obj.t(1));
            zeroPaddedTime = min(obj.t(:)) + dwell_time*((1:obj.zeroPadSize)-1)';
            individualSignals = obj.calcComponentTimeDomainSignal(zeroPaddedTime);
            fittedSignal = obj.calcTimeDomainSignal(zeroPaddedTime);
            residualSignal = obj.timeDomainSignal - fittedSignal;
            
            % Calculate lorentzian curves for each component
            nComponents = length(obj.area);
            fMat = repmat(zeroPaddedTime,[1 nComponents]);
            
            legendStrings = cell(1, nComponents);
            for iComp=1:nComponents
                legendStrings{iComp} = ['C\_' sprintf('%03.0f',iComp)];
            end
            
            % Show results to user
            ax2 = subplot(4,1,2);
            plot(zeroPaddedTime,abs(obj.timeDomainSignal),'-b');
            hold on;
            plot(zeroPaddedTime,abs(fittedSignal),'-g');
            plot(zeroPaddedTime,abs(residualSignal),'-r');
            hold off;
            ylabel('Magnitude Intensity');
            set(ax2,'xticklabel',{[]}) ;
            
%             ax3 = subplot(5,1,3);
%             phaseSig = angle(zeroPaddedTimeimeDomainSignal);
%             phaseFit = angle(fittedSignal);
%             phaseDelta = phaseSig - phaseFit;
%             plot(zeroPaddedTime, phaseSig,'-b');
%             hold on;
%             plot(zeroPaddedTime,phaseFit,'-g');
%             plot(zeroPaddedTime,phaseDelta,'-r');
%             
%             hold off;
%             ylabel('Phase (Radians)');
%             set(ax3,'xticklabel',{[]}) ;
            
            ax4 = subplot(4,1,3);
            plot(zeroPaddedTime,real(obj.timeDomainSignal),'-b');
            hold on;
            plot(zeroPaddedTime,real(fittedSignal),'-g');
            plot(zeroPaddedTime,real(residualSignal),'-r');
            hold off;
            ylabel('Real Intensity');
            set(ax4,'xticklabel',{[]});
            
            ax5 = subplot(4,1,4);
            plot(zeroPaddedTime,imag(obj.timeDomainSignal),'-b');
            hold on;
            plot(zeroPaddedTime,imag(fittedSignal),'-g');
            plot(zeroPaddedTime,imag(residualSignal),'-r');
            hold off;
            xlabel('Time');
            ylabel('Imaginary Intensity');
            legend('Measured','Fitted','Residual');
            
            %             if(~isempty(obj.fref))
            %                 % Add PPM axis
            %                 ax1ppm = subplot(5,1,1);
            %                 set(ax1ppm,'units','normalized',...
            %                     'XAxisLocation','top','YAxisLocation','right',...
            %                     'YTick',[],'YTickLabel',[],'Color','none');
            %
            %                 ax1 = axes('Position',get(ax1ppm,'Position'));
            %             else
            ax1 = subplot(4,1,1);
            %             end
            
            plot(fMat,real(individualSignals));
            legend(legendStrings);
            ylabel('Component Intensity');
            set(ax1,'xticklabel',{[]}) ;
            %             if(~isempty(obj.fref))
            %                 set(ax1ppm,'XDir','reverse');
            %             end
            
            % Keep all x axes in sinc
            linkaxes([ax1,ax2,ax4 ax5],'x');
            
            %             if(~isempty(obj.fref))
            %                 % Initialize XLim in correct units
            %                 set(ax1ppm,'xlim',NMR_Mix.relFreqToPpm(get(ax2,'XLim'),obj.fref));
            %
            %                 % Keep ppm axiz in sinc with freq axis
            %                 xLimListener = addlistener( ax1, 'XLim', 'PostSet', ...
            %                     @(src,evt) set(ax1ppm,'XLim',...
            %                     NMR_Mix.relFreqToPpm(get(ax1,'XLim'),obj.fref)) );
            %             end
        end
    end
end
