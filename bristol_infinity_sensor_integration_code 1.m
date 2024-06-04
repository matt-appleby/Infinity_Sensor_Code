function [timeVector,signalOut]=bristol_infinity_sensor_integration_code(signalIn,T_PERIOD,GAIN,GRADIENT,OFFSET)%--- is it worth calling it infinity sensor post processing as the code also applies an offset
% % Post processing function for Bristol Infinity Sensor
%  
%   Description:
%       This code is used to convert the raw signal output of a Bristol infinity sensor captured by an oscilloscope
%       from a measurement of di/dt into a current measurement. This is accomplished by applying an offset 
%       to the signal followed by numerical integration. The offset is required to compensate for the sampling error of the scope.
%       
%
%   Input variable details :
%
%       signalIn:       The double to be integrated. If this is the only
%                       input, it will be integrated with the default
%                       setting (see other variables)
%
%       T_PERIOD:       Optional argument (Integer). The sampled time period in seconds.The
%                       default T_PERIOD is 1.6*10^-10 (i.e. 0.16 ns) -- is this definetly an integer as it has decimal places??
%
%       GAIN:           Optional argument (Integer). The integration constant in V/(A/ns). The
%                       default GAIN is 1 (i.e. 1 V/(A/ns))
%
%       GRADIENT:       Optional argument (Integer). The gradient of signalOut before event in---------I would be tempted to remove GRADIENT to simplify
%                       A/us. The default GRADIENT is 0.
%
%       OFFSET:         Optional argument (Integer). The amplitude of the first value of    ----------Is it worth naming this DC_OFFSET or something else to avoid confusion with meanOffset
%                       signalOut in A The default OFFSET is 0.
%
%  Output variable details :
%
%       timeVector:     Double. The time vector of the signalOut in ns.
%
%       signalOut:      Double. The integrated signal.
%
%==========================================================================================================================================================
%
% Version history :
%
%   1.0         -- First release.
%
%============================================================================================================================================================
%
%  Authorship & Acknowledgements :
%
%   Written by Yushi Wang, Electrical Energy Management Group, University
%   of Bristol, UK. Developed with MATLAB R2021b. The code herein may
%   be freely shared and modified as long as this comment box is included
%   in its entirety and any changes and respective authors are summarised
%   here.
%
%============================================================================================================================================================

    % Initialise optional variables
    if ~exist('GAIN','var')
        GAIN=1;
    end
        GAIN=GAIN*10^9;
    if ~exist('T_PERIOD','var')
        T_PERIOD=1.6*10^-10;
    end
    if ~exist('GRADIENT','var')
        GRADIENT=0;
    end
    if ~exist('OFFSET','var')
        OFFSET=0;
    end

    % Determine the required offset by finding the average signal value where there is no event
    % 1. Determine where there is no event
    [peak,peakIndex] = max(abs(signalIn));
    signalBeforePeak=signalIn(1:peakIndex);
    for lv=1:peakIndex
        if abs(signalBeforePeak(lv))<(peak*0.05)
            signalBeforePeak(lv)=0;
        else
            signalBeforePeak(lv)=signalBeforePeak(lv);
        end
    end
    roughStartIndex=find(signalBeforePeak~=0,1);
    differentiatedSignal=flip(signalIn(1:roughStartIndex));
    startIndex=roughStartIndex-find(sign(differentiatedSignal)~=sign(differentiatedSignal(1)),1);
    % 2. Determine the average value in this period 
    meanOffset=mean(signalIn(1:round(startIndex)));

    % Create output vector
    signalLength=length(signalIn);
    timeVector=(1:signalLength)*T_PERIOD;

    % Remove the oscilloscopes sampling offset from the measurement and integrate the signal to give current.
    signalOut=cumtrapz(T_PERIOD,signalIn-meanOffset)*GAIN+GRADIENT*timeVector'*10^6+OFFSET;%%%%% Matt "What is the commer for?"
    timeVector=timeVector*10^9;
end