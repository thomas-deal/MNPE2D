function surface = RoughSurface1D(SurfType,L,NFFT,H,tide,varargin)
%% function surface = RoughSurface1D(SurfType,L,NFFT,H,tide,varargin)
% function surface = RoughSurface1D(SurfType=-1,L,NFFT,H,tide,wvheight,wvlength)
% function surface = RoughSurface1D(SurfType=0,L,NFFT,H,tide)
% function surface = RoughSurface1D(SurfType=1,L,NFFT,H,tide,seed,t,wvdir,rmsrough,corrl)
% function surface = RoughSurface1D(SurfType=2,L,NFFT,H,tide,seed,t,wvdir,wspd)
% function surface = RoughSurface1D(SurfType=3,L,NFFT,H,tide,seed,t,wvdir,wspd,fetch)
% function surface = RoughSurface1D(SurfType=4,L,NFFT,H,tide,seed,t,wvdir,f,Sf)
%
% Generates a realization of a 1D rough surface based on user inputs.
% Returns a structure with a griddedInterpolant object for evaluating the 
% surface at any point r, surface displacement eta(r), and the first two 
% derivatives of eta with respect to r.
%
% SurfType: 
%          -1 - Suinusoidal surface
%           0 - Flat surface
%           1 - RMS roughness and correlation length
%           2 - Pierson-Moskowitz fully developed seas
%           3 - JONSWAP limited fetch
%           4 - User-defined 1D spectrum
%
% Required Inputs:
%           L           - Length of surface, m
%           NFFT        - FFT size (output dimensions are NFFTx1)
%           H           - Water depth at r=0, m
%           tide        - Tide height (positive up) referenced to MLLW, m
%
% SurfType Dependent Inputs:
%           seed        - Integer seed value for repeatable random numbers
%           t           - Time since initial wave realization for 
%                         propagating waves forward in time, s
%           wvdir       - Wave direction: +1 = from +r, -1 = from -r,
%                         0 = standing waves
%           rmsrough    - RMS roughness for fully random surface, m
%           corrl       - Correlation leangth for fully random surface, m
%           wspd        - Wind speed for P-M and JONSWAP spectra, m/s
%           fetch       - Wind fetch for JONSWAP spectra, m
%           f           - Frequency for user-defined 1D spectrum, Hz
%           Sf          - User-defined spectral density, m^2/Hz
%           wvheight    - Wave height for sinusoidal surface, m
%           wvlength    - Wavelength for sinusoidal surface, m
%
% Outputs:
%           surface     - Structure with the following elements
%             .spec     - Matrix containing column of positive wavenumbers
%                         and column of one-sied spectrum magnitude
%             .f_r      - Matlab griddedInterpolant object that can be used
%                         to evaluate the surface at points r not in the
%                         original computational grid. By default it
%                         returns with the surface elevation, eta, but it
%                         can also interpolate the other derivatives. See
%                         help on griddedInterpolant.
%             .eta      - Surface height, m
%             .detadr   - First derivative of eta with respect to r, m/m
%             .d2etadr2 - Second derivative with respect to r, 1/m
%

%% Constants
g = 9.80665;
%% Vectors
dr = L/NFFT;                        % Spatial resolution, m
dk = 2*pi/L;                        % Fundamental spatial frequency, 1/m
r = (0:NFFT-1)'*dr;                 % Spatial vector, m
kx = [0:NFFT/2-1 -NFFT/2:-1]'*dk;   % Wavenumber vector, 1/m
kr = abs(kx);
kp = (0:NFFT/2-1)'*dk;              % One-sided wavenumber vector, 1/m
omega = sqrt(g*kp.*tanh(kp*H));
domegadk = min((g*tanh(kp*H)+g*kp*H.*sech(kp*H).^2)./(2*sqrt(g*kp.*tanh(kp*H))),1e3);
fk = 1/(2*pi)*sqrt(g*kp.*tanh(kp*H));
dfdk = 1/(2*pi)*(g*tanh(kp*H) + g*H*kp.*sech(kp*H).^2)./(2*sqrt(g*kp.*tanh(kp*H)));
%% Compute Spectrum
if (SurfType==-1)	% Sinusoidal Surface
    wvheight = varargin{1};
    wvlength = varargin{2};
    ksin = 2*pi/wvlength;
    eta = -tide+wvheight*sin(ksin*r);
    detadr = ksin*wvheight*cos(ksin*r);
    d2etadr2 = -ksin^2*wvheight*sin(ksin*r);
    w2 = fft(eta)/NFFT;
    w1 = 2*abs(w2(1:NFFT/2));   % Note this makes a single-sided spectrum 
    w1(1) = w2(1);              % with the same frequency and amplitude as
                                % the sinusoidal surface but the phase
                                % information is lost.                          
elseif (SurfType==0)	% Flat Surface
    eta = -tide*ones(size(r));
    detadr = zeros(size(r));
    d2etadr2 = zeros(size(r));
    w1 = zeros(size(kp));
    w1(1) = -tide;
else
    seed = varargin{1};
    t = varargin{2};
    wvdir = sign(varargin{3});
    switch SurfType
        case 1	% RMS roughnesws and correlation length
            rmsrough = varargin{4};
            corrl = varargin{5};
            w1 = (1+corrl^2*kp.^2).^(-rmsrough/2+0.5);
        case 2  % Pierson-Moskowitz Spectrum
            wspd = varargin{4};
            alpha = 0.0081;
            beta = 0.74;
            omega0 = g/wspd;
            w1 = alpha*g^2./omega.^5.*exp(-beta*(omega0./omega).^4);
            w1 = w1.*domegadk;
        case 3  % JONSWAP Spectrum
            wspd = varargin{4};
            fetch = varargin{5};
            alpha=0.076*(g*fetch/wspd)^-0.22;
            omega0=7*pi*(g/wspd)*(g*fetch/wspd^2)^-0.33;
            gamma0=3.3;
            sigma = 0.07*ones(size(omega));
            sigma(omega>omega0) = 0.09;
            delta = exp(-(omega-omega0).^2./(2*sigma.^2*omega0^2));
            w1 = alpha*g^2*gamma0.^delta./omega.^5.*exp(-1.25*(omega/omega0).^-4);
            w1 = w1.*domegadk;
        case 4  % User-Defined Spectrum
            f = varargin{4};
            Sf = varargin{5};
            Sfk = interp1(f,Sf,fk,'linear',0);
            w1 = Sfk.*dfdk;
        otherwise
            w1 = zeros(size(kp));
    end
    w1(1) = 0;
    % Two-Sided Spectrum
    switch wvdir
        case -1 % Coming from -r
            S2 = [zeros(NFFT/2+1,1); flipud(w1(2:end))]*dk;
        case 0  % Standing waves
            S2 = [w1; 0; flipud(w1(2:end))]/2*dk;
        case 1  % Coming from +r
            S2 = [w1; zeros(NFFT/2,1)]*dk;
    end
    % Random Fourier Amplitudes
    rng(seed);
    ranr = randn(NFFT,1);
    rani = randn(NFFT,1);
    ranph = 1/sqrt(2)*(ranr + 1i*rani);
    % Dispersion
    omega2 = sqrt(g*kr.*tanh(kr*H));
    % Nonzero frequencies
    eta0 = ranph.*sqrt(S2);
    etahat = zeros(NFFT,1);
    etahat(2:NFFT) = 1/sqrt(2)*(eta0(2:NFFT).*exp(-1i*omega2(2:NFFT)*t) + ...
                     conj(flipud(eta0(2:NFFT))).*exp(1i*omega2(2:NFFT)*t));
    % Compute first derivative
    etahatd = 1i*kx.*etahat;
    % Modify Nyquist value
    etahatd(NFFT/2+1) = real(etahatd(NFFT/2+1));
    % Compute second derivative
    etahatdd = -kx.^2.*etahat;
    % Convert to spatial domain
    eta = fft(etahat,NFFT);
    detadr = fft(etahatd,NFFT);
    d2etadr2 = fft(etahatdd,NFFT);
    % Subtract tide height
    eta = eta-tide;
    w1(1) = -tide;
end
%% Assign outputs
surface.spec = [kp w1];
surface.f_r = griddedInterpolant(r,real(eta),'linear','nearest');
surface.eta = real(eta);
surface.detadr = real(detadr);
surface.d2etadr2 = real(d2etadr2);
