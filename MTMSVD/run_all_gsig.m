function run_all_gsig(ncfile)

[p,~,~]=fileparts(mfilename('fullpath'));

addpath([p '/../']);

[~,outname,~]=fileparts(ncfile);

    if isempty(gcp('nocreate'));
        ncores=feature('numcores');
        parpool(ncores);
    end    
   
    stats_file=['./results/',outname,'_stats.mat'];
    load(stats_file);

    disp(['Reconstructing signals for: ' ncfile]);
    [tstats]=get_signal(stats);
    
    save(strcat('./results/',outname,'_tstats.mat'),'tstats');
    disp(['saved reconstruction results as ' outname '_tstats.mat']);

end
