function [SNR, h_rss, h_peak, Fchar, bw, Tchar, dur, Fpeak] = xoptimalsnr(h,t0,fs,S,F0,dF,Fmin,Fmax)
% XOPTIMALSNR - Compute SNR and other properties of a signal in noise.
%
% XOPTIMALSNR - Compute the SNR and time-frequency measures of a waveform
% in a specified noise background.  Simpler noise-independent measures of 
% the wave amplitude are also provided.
%
%   [SNR, h_rss, h_peak, Fchar, bw, Tchar, dur] = ... 
%       xoptimalsnr(h,t0,fs,S,F0,dF,Fmin,Fmax)
%
%   h      Array.  Waveform timeseries data.  Each column holds the 
%          timeseries for one of the GW polarizations.  (There may be any
%          number of polarisations, for the non-GR case.  The order of 
%          polarizations does not matter.)
%   t0     Scalar.  Time at which the first waveform data point h(1,:) is
%          sampled.   
%   fs     Scalar.  Sampling rate (Hz) of waveform data.
%   S      Vector (optional).  Noise background one-sided POWER (not
%          amplitude) spectrum.  If supplied, it must cover at least the
%          range [Fmin, Fmax], and must be linearly sampled in ascending
%          order of frequency.  For noise-independent signal measures use
%          S=f0=df=[].  In this case the constant spectrum S(f)=2 will be
%          used (corresponding to a two-sided noise spectrum of unity).
%   F0     Scalar (optional).  Frequency at which S(1) is sampled.
%   dF     Scalar (optional).  Frequency spacing of the noise spectrum S.
%   Fmin   Scalar (optional).  Minimum frequency (Hz) to include in 
%          frequency-domain calculations.
%   Fmax   Scalar (optional).  Maximum frequency (Hz) to include in 
%          frequency-domain calculations.
%
% Computations are done in the time domain (TD) and frequency
% domain (FD) using the energy distributions
%      p_TD = h(:,1).^2 + h(:,2).^2;
%      p_FD = 2(|FFT(h(:,1))|.^2 + |FFT(h(:,2))|.^2);  % for f>=0
% With these conventions, the output is
%   SNR    Signal-to-noise ratio of h in the given noise background, 
%          defined as 
%            SNR = (2 \int_Fmin^Fmax df p_FD./S).^0.5
%   h_rss  The root-sum-square amplitude (Hz^-0.5) of the waveform h:
%            h_rss = \int_Fmin^Fmax df p_FD
%   h_peak The maximum absolute amplitude of the waveform h:
%            h_peak = max(p_TD).^0.5
%   Fchar  Characteristic frequency (Hz) of h in the given noise  
%          background, defined as 
%                    \int_Fmin^Fmax df f p_FD./S
%            Fchar = ----------------------------------------
%                     \int_Fmin^Fmax df p_FD./S
%          where Fmin = max(f) and \tilde(h) is the FFT of h.  
%   bw     Effective bandwidth (Hz) of h in the given noise  
%          background, defined as 
%                 \int_Fmin^Fmax df (f-Fchar).^2 p_FD./S
%            bw = ---------------------------------------------------
%                  \int_Fmin^Fmax df p_FD./S
%          where Fmin = max(f) and \tilde(h) is the FFT of h.  
%   Tchar  Characteristic time at which the event occurs, defined as 
%                    \int dt t p_TD(t)
%            Tchar = -----------------
%                     \int dt p_TD(t)
%   dur    Duration of the signal, defined as 
%                    \int dt (t-Tchar).^2 p_TD(t)
%            Tchar = ----------------------------
%                     \int dt p_TD(t)
%
% Note that no windowing is used for computing FFTs; windowing and 
% averaging is not really sensible for a transient signal, since by
% definition the statistical properties of the signal are changing over
% time.  There may be problems for, e.g., band-passed noise bursts.
%
% Note that the SNR and h_rss measures are restricted to the frequency
% interval [Fmin,Fmax], h_peak is evaluated in the time domain (i.e., using
% data from the full frequency range. 
%
% Do "type FourierTransformConventions" for information on the conventions
% used for Fourier transforms.
% 
% $Id: xoptimalsnr.m 3911 2012-07-03 13:51:31Z patrick.sutton@LIGO.ORG $


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Checks.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ---- Valid number of input arguments.
error(nargchk(3,8,nargin,'struct'))

% ---- Optional arguments.
if ( (nargin<8) || isequal(Fmax,[]) )
    Fmax = floor(fs/2);
end
if ( (nargin<7) || isequal(Fmin,[]) )
    Fmin = 0;
end
% ---- Is there an input noise spectrum?
if ( (nargin<6) || isequal(S,[]) || isequal(F0,[]) || isequal(dF,[]) )
    % ---- Set flag to make dummy noise spectrum.
    noise = 0;
else
    noise = 1;
    % ---- Validate noise spectrum data.
    if (~isscalar(dF) | dF<=0)
        error('Sprectrum frequency resolution dF must be a positive scalar.');
    end
    if (~isscalar(F0) | F0<0)
        error('Sprectrum lowest frequency F0 must be a non-negative scalar.');
    end
    if ~isvector(S)
        error('Spectrum S be a vector or empty array.');
    end
end
if (~isscalar(Fmax) | Fmax<=0)
    error('Frequency limit Fmax must be a positive scalar.');
end
if (~isscalar(Fmin) | Fmin<0)
    error('Frequency limit Fmin must be a non-negative scalar.');
end
if Fmin>=Fmax
    error('Frequency limits must satisfy Fmin<Fmax.');
end
% ---- Require positive sampling rate.
if fs<=0
    error('Sampling rate fs must be positive.');
end
if ~isscalar(t0)
    error('Timeseries start time t0 must be a scalar.');
end
% % ---- Remove this check to allow for non-GR polarisations.
% % ---- Input timeseries must be one- or two-column array.
% if (size(h,2)>2)
%     error('Input timeseries h must be a one- or two-column array.');
% end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Preparations.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ---- Number of data points.  Force to be even.
N = size(h,1);
if ~iseven(N)
    %warning(['Function not guaranteed for waveforms with odd number ' ...
    %    'of samples.  Dropping last sample']);
    h(end,:) = [];
    N = size(h,1);    
end

% ---- Duration of tseries.
T = N/fs;  

% ---- Sample times.
t = t0+1/fs*[0:(N-1)]';

% ---- If no noise spectrum is supplied then make a dummy noise vector 
%      covering [0,Nyquist] Hz.  This will allow us to assume henceforth
%      that S, F0, and dF are defined.
if (~noise)
    dF = 1/T;
    F0 = 0;
    S = 2*ones(N/2+1,1);
end

% ---- Force column vectors.
S = S(:);

% ---- Vector of sampled noise frequencies.
F = F0+[0:length(S)-1]'*dF;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Frequency-domain calculations.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ---- Signal.

% ---- Standard FFT; works by column.
hf = 1/fs*fft(h);
% ---- Vector holding frequencies in usual screwy FFT order:
%        vector element:   [ 1  2  ...  N/2-1  N/2    N/2+1            N/2+2   ... N-1  N ]
%        frequency (df):   [ 0  1  ...  N/2-2  N/2-1  (N/2 or -N/2)   -N/2+1  ... -2   -1 ]
f = [ 0:N/2 , -N/2+1:-1 ]'/T;
% ---- Keep only frequency components in [Fmin,Fmax].  Note that most of
%      the frequency-domain formulas include a factor of 2 to account for
%      negative frequencies. 
index = find(f>=Fmin & f<=Fmax);
f = f(index);
% ---- Distribution of signal energy in frequency.
p_FD = 2*sum(hf(index,:).*conj(hf(index,:)),2);
% ---- f=0 and f=Nyq bins should count for 1/2.
if (f(1)==0)
    p_FD(1) = 0.5*p_FD(1);
end
if (f(end)==fs/2)
    p_FD(end) = 0.5*p_FD(end);
end

% ---- Noise.

% ---- Does vector of sampled noise frequencies cover requested range?
if ( (F(1)>f(1)) | (F(end)<f(end)) ) 
    error('Noise spectrum does not cover desired frequency range.')
end
% ---- Force interpolation of S from sampled frequencies F to data
%      frequencies f.
S = exp(interp1(log(F),log(S),log(f)));
F = f;
dF = 1/T;
F0 = F(1);


% ---- All set to do calculations.  Assured that f, F interchangable.

% ---- SNR^2 versus frequency.
SNR2f = 2*p_FD./S;

% ---- SNR on [Fmin,Fmax].
SNR = (dF*sum(SNR2f)).^0.5;

% ---- Characteristic frequency.
Fchar = sum(f.*SNR2f)./sum(SNR2f);

% ---- Peak Frequency
Fpeak = f(p_FD==max(p_FD));
if length(Fpeak)>1
	Fpeak=Fpeak(1);
end

% ---- Characteristic bandwidth.  
bw = (sum((f-Fchar).^2.*SNR2f)./sum(SNR2f))^0.5;

% % ------ Frequency interval containing central 50% of total |h(f)|^2.
% % ------ Accumulation of |h(f)|^2 over frequency.
% chf2 = cumsum(SNR2f)./sum(SNR2f);
% % ------ Select out unique points (necessary for linear interp)  
% [Z,I] = unique(chf2);
% chf2 = chf2(I);
% zf = f(I);
% f75 = interp1(chf2,zf,0.75);
% f25 = interp1(chf2,zf,0.25);
% bw = f75 - f25;

% ---- RSS amplitude.
h_rss = (dF*sum(p_FD)).^0.5;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Time-domain calculations.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ---- Distribution of signal energy in time.
p_TD = sum(h.^2,2);

% ---- Peak amplitude.
h_peak = max(p_TD).^0.5;

% ---- Characteristic time.
Tchar = sum(t.*p_TD)./sum(p_TD);

% ---- Duration.  
dur = (sum((t-Tchar).^2.*p_TD)./sum(p_TD))^0.5;

% % ------ Time containing 50% of total h(t)^2.
% % ------ Accumulation of h(t)^2 over time.
% cht2 = cumsum(h.^2)/(h'*h);
% % ------ Select out unique points (necessary for linear interp)  
% [Z,I] = UNIQUE(cht2);
% cht2 = cht2(I);
% zt = t(I);
% te = interp1(cht2,zt,0.75);
% ts = interp1(cht2,zt,0.25);
% dur = te - ts;


% ---- Done
return
