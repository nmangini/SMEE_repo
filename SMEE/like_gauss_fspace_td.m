function logL = like_gauss_fspace(wave_ft, noise, deltaF, len, freqs, T_shift, model,varargin)
% This function calculates the likelihood, assuming a Gaussian
% distribution, given some data (containing a signal and noise), the
% PSD of the noise, and a model function.
%
% The code will return the likelihood value, and the natural logarithm of
% the likelihood value.

% check whether model is a function handle (or string that is a function
% name)
if ischar(model)
    f = str2func(model);
elseif isa(model, 'function_handle')
    f = model;
else
    error('Error... Expecting a model function!');
end

% evaluate the model, including a time shift
md=feval(f, varargin{:});
md_ft = md;%.*exp(-2*pi*1i*freqs*T_shift);

% check that the length of the data and the model are the same
if length(wave_ft) ~= length(md_ft)
    error('Error... Length of data and model are not the same!');
end

% get the log likelihood
%logL = -2*(1/(deltaT*len))*sum(((abs(wave_ft - md_ft)).^2)./(abs(noise).^2)); ...
logL = -2*deltaF*sum(((abs(wave_ft(30:end) - md_ft(30:end))).^2)./(noise(30:end)));
%logL = -sum(((abs(wave_ft - md_ft)).^2)./(abs(noise).^2));
%logL = -sum(((abs(wave_ft - md_ft)).^2)./noise_PSD); ...
    %- 0.5*sum(log(2*deltaT/(pi*len*Pf1)));

% get likelihood
L = exp(logL);
