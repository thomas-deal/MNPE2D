function alphaout = ConvertAttenuationUnits(alphain,fromunits,tounits,varargin)
%% function alphaout = ConvertAttenuationUnits(alphain,fromunits,tounits)
% function alphaout = ConvertAttenuationUnits(alphain,fromunits,tounits,f,c)
%
% Converts attenuation from one set of units to another. Units are 1/m,
% dB/km, dB/m/kHz, dB/km/Hz, dB/lambda
%
% Inputs:
%       alphain     - Input attenuation
%       fromunits   - Units of alphain
%       tounits     - Units of alphaout
%
% Optional Inputs:
%       f           - Frequency, Hz
%       c           - Sound speed, m/s
%
% Outputs:
%       alphaout    - Output attenuation
%

%% Check Input Arguments
% Default values
deff = true;
defc = true;
f = 1e3;
c = 1500;
switch nargin
    case 4
        deff = false;
        f = varargin{1};
    case 5
        deff = false;
        defc = false;
        f = varargin{1};
        c = varargin{2}; 
end
e = exp(1);
%% Convert to 1/m
switch fromunits
    case '1/m'
        alpha = alphain;
    case 'dB/km'
        alpha = alphain/(20e3*log10(e));
    case 'dB/m/kHz'
        if deff
            disp(['Frequency not specified - using default value of ' num2str(f) ' Hz']);
        end
        alpha = alphain*f/(20e3*log10(e));
    case 'dB/km/Hz'
        if deff
            disp(['Frequency not specified - using default value of ' num2str(f) ' Hz']);
        end
        alpha = alphain*f/(20e3*log10(e));
    case 'dB/lambda'
        if deff && defc
            disp(['Frequency and sound speed not specified - using default values of ' num2str(f) ' Hz and ' num2str(c) ' m/s']);
        elseif defc
            disp(['Sound speed not specified - using default value of ' num2str(c) ' m/s']);
        end
        alpha = alphain*f/(20*log10(e)*c);
    otherwise
        disp(['Error: input units ' fromunits ' not recognized. Returning input alpha.'])
        alphaout = alphain;
        return
end
%% Convert to Output Units
switch tounits
    case '1/m'
        alphaout = alpha;
    case 'dB/km'
        alphaout = 20e3*log10(e)*alpha;
    case 'dB/m/kHz'
        if deff
            disp(['Frequency not specified - using default value of ' num2str(f) ' Hz']);
        end
        alphaout = 20e3*log10(e)/f*alpha;
    case 'dB/km/Hz'
        if deff
            disp(['Frequency not specified - using default value of ' num2str(f) ' Hz']);
        end
        alphaout = 20e3*log10(e)/f*alpha;
    case 'dB/lambda'
        if deff && defc
            disp(['Frequency and sound speed not specified - using default values of ' num2str(f) ' Hz and ' num2str(c) ' m/s']);
        elseif defc
            disp(['Sound speed not specified - using default value of ' num2str(c) ' m/s']);
        end
        alphaout = 20*log10(e)*c/f*alpha;
    otherwise
        disp(['Error: output units ' tounits ' not recognized. Returning units of 1/m.'])
        alphaout = alpha;
end
        