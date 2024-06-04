function [timeVector,signalOut]=bristol_infinity_sensor_integration_code(signalIn,T_PERIOD,GAIN,GRADIENT,OFFSET)
%   [timeVector,signalOut]=bristol_infinity_sensor_integration_code(signalIn,[T_PERIOD],[GAIN],[GRADIENT_CORRECTION],[OFFSET_CORRECTION])
%
%   Where input variables in square brackets are optional - i.e. The
%   function accepts anywhere from 1 to 5 inputs, inclusive. This code
%   integrates the input signal.
%
%% ==========================================================================================================================================================
%
%% Version history :
%
%   1.0         -- First release.
%
%============================================================================================================================================================
%
%% Notes :
%
%   The code uses the following naming convention :
%       _____________________________________________
%                         |
%           Item(s)       |   Naming convention
%       __________________|__________________________
%           Variables     |   initialLowerTitleCase
%       ------------------+--------------------------
%           Functions     |   lower_case(...)
%       ------------------+--------------------------
%           constants     |   UPPER_CASE
%       __________________|__________________________
%
%============================================================================================================================================================
%
%% Input variable details :
%
%   Input variables: 
%
%       signalIn:       The double to be integrated. If this is the only
%                       input, it will be integrated with the default
%                       setting (see other variables)
%
%       GAIN:           Integer. The integration constant in V/(A/ns). The
%                       default GAIN is 1 (i.e. 1 V/(A/ns))
%
%       T_PERIOD:       Integer. The sampled time period in seconds.The
%                       default T_PERIOD is 1.6*10^-10 (i.e. 0.16 ns)
%
%       GRADIENT:       Integer. The gradient of signalOut before event in
%                       A/us. The default GRADIENT is 0.
%
%       OFFSET:         Integer. The amplitude of the first value of
%                       signalOut in A The default OFFSET is 0.
%% Output variable details :
%
%       timeVector:     Double. The time vector of the signalOut in ns.
%
%       signalOut:      Double. The integrated signal.
%
%============================================================================================================================================================
%
%% Authorship & Acknowledgements :
%
%   Written by Yushi Wang, Electrical Energy Management Group, University
%   of Bristol, UK. Developed with MATLAB R2021b running on Windows 11. The
%   structure of the code is taken from Harry Dymond. The code herein may
%   be freely shared and modified as long as this comment box is included
%   in its entirety and any changes and respective authors are summarised
%   here.
%
%% ==========================================================================================================================================================

%% Initialise optional vartiable 
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

%% Find event rising edge index
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

%% Create output vector
    signalLength=length(signalIn);
    timeVector=(1:signalLength)*T_PERIOD;
    meanOffset=mean(signalIn(1:round(startIndex)));
    signalOut=cumtrapz(T_PERIOD,signalIn-meanOffset)*GAIN+GRADIENT*timeVector'*10^6+OFFSET;
    timeVector=timeVector*10^9;
    % plot(timeVector,signalOut)
    % hold on
end