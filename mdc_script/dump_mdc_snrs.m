% dump_mdc_snrs.m
%
% Use xoptimalsnr.m to produce some useful characterisation data from an ascii
% dump of an MDC frame.
%
% USAGE:
%
%	dump_mdc_snrs(injections,times,ifo_site,era)
%
%	INPUT
%		injectionfile = NINJA MDC frame
%		channel = channel to extract from frame
%		times = ascii dump of times from sim_inspiral table (expects 4 columns:
%		geocentric, h_end_time, l_end_time, v_end_time)
%		ifo_site = detector location (one of 'h','l' or 'v')
%		era = 'early' or 'nominal' noise curve
%
% 
% James Clark, July 2012

function dump_mdc_snrs(injectionfile,channel,timesfile,Fmin,Fmax,era)

% --- Paths & hardcoded variables
SMEEBBH_PREFIX=getenv('SMEEBBH_PREFIX');
srdpath=(sprintf('%s/SRDs/',SMEEBBH_PREFIX));
fs=16384;
injscale=1;
dF=0.5;

% parse filename
[sites,tag,gpsStart,duration]=strread(injectionfile,...
	'%s %s %f %f','delimiter','-');
[ifo,type]=strread(channel,'%s %s','delimiter',':');
ifo=ifo{1};

% --- Logic to parse input
nominal_aligo='T0900288_V3_ZERO_DET_high_P.txt';
%early_aligo='early_aligo.dat';
early_aligo='T1200307_V4_early_aligo_asd.txt';
AdvVirgo='design_virgo_T1300121_V1.txt';

% default era is nominal advanced era
if nargin<6 | isempty(era)
	era='nominal';
end

% check validity of era
switch era
case 'nominal'
	disp('Using nominal design (~2020) sensitivity noise curve')
case 'early'
	disp('Using estimated early (~2015) sensitivity noise curve')
otherwise
	error('invalid era - should be one of "nominal" or "early" or unspecified')
end

% read times file - reading strings is an easy way to preserve precision
disp('reading times file')
[geocentric_times,h_times,l_times,v_times,...
	distance, eff_dist_h, eff_dist_l, eff_dist_v]=textread(timesfile,...
		'%s %s %s %s %s %s %s %s','delimiter',',');

% decide which noise curve to load
switch ifo

case {'H1','L1'}
	switch era
	case 'nominal'
		noise=load([srdpath,nominal_aligo]);
		noise_freqs=noise(:,1);
		noise_asd=noise(:,2);
	case 'early'
		noise=load([srdpath,early_aligo]);
		noise_freqs=noise(:,1);
		noise_asd=noise(:,2);
	end

	switch ifo
	case 'H1'
		injtimes=h_times;
		effdist=eff_dist_h;
	case 'L1'
		injtimes=l_times;
		effdist=eff_dist_l;
	end

case 'V1'
		noise=load([srdpath,AdvVirgo]);
		noise_freqs=noise(:,1);
	switch era
	case 'nominal'
		noise_asd=noise(:,2);
	case 'early'
		noise_asd=noise(:,2);
	end
	injtimes=v_times;
	effdist=eff_dist_v;

otherwise
	error('invalid detector site - should be one of "h", "l" or "v"')
end
	
% read injection file
disp('reading injection file...')
[signal_hoft,signal_times]=frgetvect(injectionfile,channel,gpsStart,duration);
signal_times = signal_times + gpsStart;

disp('injection data read')

% --- Interpolate the noise spectrum so we have uniform frequency sampling
noise_freq_interp=10:dF:8192;
noise_asd_interp=exp(interp1(log(noise_freqs), log(noise_asd), ...
	log(noise_freq_interp)));

% --- Computations

% reduce injection list to the current frame
injtimes_num=str2num(cell2mat(injtimes));
injtimes=injtimes(injtimes_num>=gpsStart & injtimes_num <gpsStart+duration);

% open results file
outfilename=strrep(injectionfile,'.gwf',sprintf('-%s_characteristics.dat',era));
outfilename=strrep(outfilename,'GHLTV',ifo);
fid=fopen(outfilename,'w');

% Loop through injection times; if there exist signal timestamps within 200 ms
% of each injection time, pass the data to xoptimalsnr().  This condition lets
% us use timestamp files which contain a superset of injection times to those
% actually in the data file being used

disp('beginning loop over injections')
for i=1:length(injtimes)

	if sum(signal_times>=str2double(injtimes{i})-0.5 & signal_times<str2double(injtimes{i})+0.5)

	   injection=signal_hoft( signal_times>=str2double(injtimes{i})-0.5 ...
			& signal_times<str2double(injtimes{i})+0.5 );

	   [snr,hrss,hpeak,fchar,bw,tchar,dur,fpeak]=xoptimalsnr(injection, ...
			0, fs, noise_asd_interp.^2 ,10 , dF, Fmin, Fmax);

	   disp(sprintf('%d/%d: GPSstart=%s, phys distance=%s, eff dist=%s, snr=%1.2f, hrss=%1.2e, fchar=%1.2f, bw=%1.2f, dur=%1.2f (ms), fpeak=%1.2f', ...
		   i, length(injtimes), injtimes{i}, distance{i}, effdist{i}, snr, hrss, fchar, bw, 1000*dur, fpeak))

		fprintf(fid,'%s %s %s %.10f %.10e %.10f %.10f %.10f %.10f\n', ...
			injtimes{i}, distance{i}, effdist{i}, snr, hrss, fchar, bw, dur, fpeak);
	end

end % loop over injections
fclose(fid)

disp('...finished!')

exit

