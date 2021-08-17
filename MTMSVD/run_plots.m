function run_plots(ncfile,lf,rf)

    [~,outname,~]=fileparts(ncfile);
    
    [p,~,~]=fileparts(mfilename('fullpath'));

    addpath([p '/../']);

    stat_file=[p '/results/' outname '_tstats.mat'];
    clv_file=[p '/results/' outname '_clv.mat'];

    tstats=load(stat_file);
    tstats=tstats.tstats;
    load(clv_file);
    
    itrdb=nc2struct(ncfile);
    [~,bs]=trim_chrons(itrdb);
    
    ncdat.lat=itrdb.lat.data;
    ncdat.lat(bs)=[];
    ncdat.lon=itrdb.lon.data;
    ncdat.lon(bs)=[];
    
    %disp('Making confidence plot');
    %conf_name=strcat(p,'/results/',outname,'_conf.png');
    %confonly(tstats,clvstat.lfvall,100,clvstat.med_lfv,clvstat.med_f0,clvstat.lfv_set,conf_name);
   
    lfs=num2str(lf);
    lfs=lfs(3:end);
    rfs=num2str(rf);
    rfs=rfs(3:end);

    map_name=strcat(p,'/results/',outname,'_',lfs,'-',rfs,'_map.png');
    disp(['Making quiver plot: ',num2str(lf),'-',num2str(rf)]);
    
    [~,f0,~]=get_peaks(tstats,lf,rf);
    mkMTMSVDplots(ncfile,map_name,ncdat,f0);

    disp(['Making time series plot: ',num2str(lf),'-',num2str(rf)]);
    ts_plot(tstats,outname,lf,rf);

end    
