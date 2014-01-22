function load_frame(catalogue)
% This is a test of how to load frame files using for-loops and the
% command frgetvect
warning('off','frgetvect:info')

tic
catalogue = 'Q';

if strcmp('Q',catalogue) == 1
    L_catalogue = 13;
    injections = zeros(49152,L_catalogue);
    qvalues = [1; 115; 130; 145; 150; 160; 175; 190; 2; 205; 220; 235; 250];
    for i = 1:L_catalogue
    fname = strcat('/data/nmangini/SMEE_BBH/Waveforms/MDCframes/',catalogue,...
        '-series/GHLTV-',catalogue,'-series_q',num2str(qvalues(i)),'-900000005-105.gwf');
    [injections(:,i),ts] = frgetvect(fname,'H1:STRAIN',900000013,3);
    end
elseif strcmp('HR',catalogue) == 1
    L_catalogue = 15;
    injections = zeros(49152,L_catalogue);
    qvalues = [1; 1; 1; 1; 1; 1; 1; 1; 1; 150; 150; 150; 2; 250; 4];
    avalues = [0; 01; 02; 03; 04; 05; 06; 07; 08; 0; 02; 04; 0; 0; 0];
    for i = 1:L_catalogue
    fname = strcat('/data/nmangini/SMEE_BBH/Waveforms/MDCframes/',catalogue,...
        '-series/GHLTV-',catalogue,'series_q',num2str(qvalues(i)),'_a',num2str(avalues(i),'%02i'),'-900000005-105.gwf');
    [injections(:,i),ts] = frgetvect(fname,'H1:STRAIN',900000013,3);
    end
elseif strcmp('RO3',catalogue) == 1
    L_catalogue = 20;
    injections = zeros(49152,L_catalogue);
    qvalues = [150; 150; 150; 150; 2; 2; 2; 2; 2; 2; 2; 250; 250; 250; 250; 250; 250; 4; 4; 4];
    avalues = [040; 060; 060; 060; 040; 060; 060; 060; 060; 060; 060; 040; 040; 040; 060; 060; 060; 060; 060; 060];
    othvalues = [060; 045; 060; 090; 060; 045; 060; 090; 135; 180; 270; 045; 060; 090; 045; 060; 090; 045; 060; 090];
    for i = 1:L_catalogue
    fname = strcat('/data/nmangini/SMEE_BBH/Waveforms/MDCframes/',catalogue,...
        '-series/GHLTV-',catalogue,'series_q',num2str(qvalues(i)),'_a',num2str(avalues(i),'%03i'),...
        '_oth',num2str(othvalues(i),'%03i'),'-900000005-105.gwf');
    [injections(:,i),ts] = frgetvect(fname,'H1:STRAIN',900000013,3);
    end
end

% Save
injection_savefile = strcat('injection_',catalogue,'-series');
save(injection_savefile,'injections')
time_savefile = strcat('Tinjection_',catalogue,'-series');
save(time_savefile,'ts');
toc
