function[snr]=computeSNR_colourednoise(sigwonoise, DetectorNoiseFile)
% Script to calculate values to characterize the intensity of the signal
% a) SNR

% Formulae:
% SNR^2 = 2 * int_{-inf}^{inf}abs(h(f))^2/S(f) df (eq. 79 Martin's notes)
%       =(4 * deltaT / N) Sum_{k} |h(f_{k})|^2/S(f_{k})
%

sig=sigwonoise;
len=length(sig);

deltaT = 1/16384;
duration=len*deltaT;
deltaF = 1/duration;
N=len;%1/(deltaT*deltaF);

if ( mod(N,2)== 0 )
  numFreqs = N/2 - 1;
else
  numFreqs = (N-1)/2;
end

f = deltaF*[1:1:numFreqs]';

sig_ft = fft(sig)*deltaT;
sig_ft = sig_ft(2:numFreqs+1);
% sig_ft = sig_ft(2:numFreqs);
%DetectorNoiseFile = DetectorNoiseFile(2:numFreqs+1);

snr = sqrt((4*(deltaF)*sum(abs(sig_ft).^2./abs(DetectorNoiseFile))));
%snr = sqrt((4*((deltaT/len))*sum(((abs(sig_ft(30:end))).^2)./((DetectorNoiseFile(30:end)).^2))));

% % noise addition
% %noise=(variance)*randn(35000,1);
% 
% sigwinoise = (sigwonoise + noise);%*1e22;
% 
% % DFT using 'pwelch' command
% nfft = fs;
% N = nfft;
% % typeofwindow = 'tukey';% r = 0.5; % wvalues = tukeywin(N,r);
% typeofwindow = 'Hamming'; % by default if window is not specified
% overlap = 0.75; noverlap = round(overlap * N);
% % calculate normalisation factors
% % S1 = sum(wvalues);
% % S2 = sum(wvalues.^2);
% % NENBW = N*S2/S1^2;
% % ENBW = NENBW * fs; % Effective Noise BandWidth
% 
% 
% [sigwiPxx,fwipwelch] = pwelch(sigwonoise,N,noverlap,nfft,fs);
% sigwiPSD = sigwiPxx; sigwiLSD = sqrt(sigwiPSD);
% 
% % Calculation of Noise Spectral Density
% [Pxxn,fn] = pwelch(noise,N,noverlap,nfft,fs);
% PSDn = Pxxn; LSDn = sqrt(PSDn);
% meanPSDn = mean(PSDn); meanLSDn = mean(LSDn);
% 
% % Calculation of SNR
% % following eq.79 from Martin's notes.
% df = 1; % frequency resolution is 1Hz because we used nfft = fs = 16384;
% Snf = LSDn;
% sighf = sigwiLSD;
% 
% snr = sqrt(4*sum((abs(sighf)).^2 ./ Snf) * df);



















