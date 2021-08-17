function run_all_svd(ncfile)

[p,~,~]=fileparts(mfilename('fullpath'));

addpath([p '/../']);

[~,outname,~]=fileparts(ncfile);

    if isempty(gcp('nocreate'));
        ncores=feature('numcores');
        parpool(ncores);
    end    
    
    itrdb=nc2struct(ncfile);

    [trim_stats,bs]=trim_chrons(itrdb);
    
    disp(['Running SVD on itrdb data: ' ncfile]);
    
    [lfv,stats]=parsvd(trim_stats.X);
    stats.lon=trim_stats.lon;
    stats.lat=trim_stats.lat;
    %save(strcat('./results/',outname,'_stats.mat'),'stats');
    %disp(['saved SVD results as ' outname '_stats.mat']);

    disp('Running analysis of itrdb SVD');
    
    conf_name=strcat('./results/',outname,'_conf.png');
    %map_name=strcat('./results/',outname,'_map.png');
    stats=parsig(stats,conf_name);
    save(strcat('./results/',outname,'_stats.mat'),'stats');
    disp(['saved SVD and analysis results as ' outname '_stats.mat']);

    %ncdat.lat=itrdb.lat.data;
    %ncdat.lat(bs)=[];
    %ncdat.lon=itrdb.lon.data;
    %ncdat.lon(bs)=[];

    %disp('Making quiver plot');
    %mkMTMSVDplots(ncfile,map_name,ncdat);

end
