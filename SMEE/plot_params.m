%plot_params Used to plot variables in SMEE output file

% Get directory names
MainDir = cd;
d = dir(pwd);
isub = [d(:).isdir];
nameFolds = {d(isub).name}';
nameFolds(ismember(nameFolds,{'.','..'})) = [];

% Loop through directories and plot parameter
for i = 1:length(nameFolds)
    cd(nameFolds{i})
    %fprintf('Looking in %s\n',nameFolds{i})
    SMEE_file = dir('*.mat');
    load(SMEE_file.name)
    plot(postbetas(:,1),'.')
    xlabel('Nested Sampling Step','fontweight','bold')
    ylabel('Beta Value','fontweight','bold')
    t = strcat(nameFolds{i},'Posterior Beta Values');
    title(t,'fontweight','bold')
    figtitle = strcat(nameFolds{i},'_postbetas');
    saveas(gcf,figtitle,'png')
    close('all')
    cd(MainDir)
end

cd(MainDir)