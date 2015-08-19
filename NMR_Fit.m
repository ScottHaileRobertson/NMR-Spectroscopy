%% NMR_Fit
% A class that fits FID to a series of exponentially decaying components using a
% NMR_Mix object. This class knows how to fit with
% or without constraints.
%
classdef NMR_Fit < NMR_Mix
    
    properties
        ub;
        lb;
        t;
        f;
        timeDomainSignal;
        spectralDomainSignal;
        zeroPadSize;
        lineBroadening;
    end
    
    methods
        function obj = NMR_Fit(time_domain_signal, t, ...
                area, freq, fwhm, phase, line_broadening, zeroPadSize)
            % Constructs a fitting object. Note - this does not perform the
            % fitting.
            
            % Construct NMR_Mix object
            obj = obj@NMR_Mix(area,freq,fwhm,phase);
            
            obj.t = t(:); % make a column vector
            nSamples = length(obj.t);
            
            % Default line broadening is zero
            obj.lineBroadening = line_broadening;
            if(isempty(obj.lineBroadening))
                obj.lineBroadening = 0;
            end
            
            % Apply linebroadening
            time_domain_signal = time_domain_signal(:); % make a column vector
            obj.timeDomainSignal = time_domain_signal.*exp(-pi*obj.lineBroadening*obj.t); % Note that pi*fwhm = 1/t2star
            
            % Default zero padding is none
            obj.zeroPadSize = zeroPadSize;
            if(isempty(obj.zeroPadSize))
                obj.zeroPadSize = nSamples;
            end
            
            % Calculate spectrum
            dwell_time = (obj.t(2)-obj.t(1));
            obj.spectralDomainSignal = dwell_time...
                *fftshift(fft(obj.timeDomainSignal,obj.zeroPadSize));
            
            % Calculate freq samples
            obj.f = linspace(-0.5,0.5,obj.zeroPadSize+1)/dwell_time;
            obj.f = obj.f(1:(end-1)); % Take off last sample to have nSamples  
            obj.f = obj.f(:); % Make a column vector
        end
        
        function obj = setBounds(obj, area_lb, area_ub, freq_lb, freq_ub,...
                fwhm_lb, fwhm_ub, phase_lb, phase_ub)
            % Sets the constraints for fitting this NMR_Fit object. 
            % Adding additional components to the NMR_Fit will erase the
            % bounds, so only add bounds immediately before you fit the 
            % signal..
            obj.ub = [area_ub(:)'; freq_ub(:)'; fwhm_ub(:)'; phase_ub(:)'];
            obj.lb = [area_lb(:)'; freq_lb(:)'; fwhm_lb(:)'; phase_lb(:)'];
            
            % Check that constraints are possible
            if(any(obj.ub < obj.lb))
                error('Impossible constraint set.');
            end
        end
        
        function describe(obj)
            % This function will report on the fitted component
            % intensities, frequencies, linewidths, and phases. Note that
            % if line broadening is used, this method subtracts off the
            % line broadening from the fitted values, thus reporting them
            % as if there was no line broadening applied
            disp('area (arbs)  Freq (Hz)  Linewidth(Hz)  Phase(degrees)');
            for iComp = 1:length(obj.area)
                disp([sprintf('%8.3e',obj.area(iComp)) ' ' ...
                    sprintf('%+8.2f',obj.freq(iComp))  '  ' ...
                    sprintf('%8.2f',obj.fwhm(iComp)-obj.lineBroadening) '  ' ...
                    sprintf('%+9.2f',obj.phase(iComp))]);
            end
        end
    end
end
