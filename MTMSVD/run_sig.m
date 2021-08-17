function run_sig(ncfile)

    [~,outname,~]=fileparts(ncfile);
    
    [p,~,~]=fileparts(mfilename('fullpath'));

    addpath([p '/../']);

    stat_file=[p '/results/' outname '_tstats.mat'];

    itrdb=nc2struct(ncfile);

    [cstats,~]=trim_chrons(itrdb);
    [~,stats]=parsvd(cstats.X);
    disp(['finished running parsvd on ', ncfile]);
    
    [tstats]=get_signal(stats);
    save(stat_file,'tstats');
    disp(['saved signal stats as',stat_file]);

end    
