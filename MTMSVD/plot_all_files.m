function plot_all_files(syear,eyear,pidx,lf,rf)

[p,~,~]=fileparts(mfilename('fullpath'));

matfile_path=[p '/results/'];

files=dir([matfile_path '/*_' num2str(syear) '-' num2str(eyear) '_stats.mat']);

for i=1:length(files)
    file=files(i).name;
    file=strcat(matfile_path,'/',file);
    dat=load(file);
    %plot_signals(dat.stats,file(1:end-10),syear,eyear,pidx,lf,rf);
    [lfv,yrs,freqs]=mtm_svd_evol_lfv(dat.stats.X,1,40);
    plot_name=[file(1:end-10) '_evol_lfv.png'];
    plot_evol_lfv(yrs,freqs,lfv,plot_name);
end    
