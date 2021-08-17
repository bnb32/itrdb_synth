function run_all_files(syear,eyear)

[p,~,~]=fileparts(mfilename('fullpath'));

ncfile_path=[p '/../results/'];

files=dir([ncfile_path '/*_' num2str(syear) '-' num2str(eyear) '.nc']);

for i=1:length(files)
    %run_all_gsig([ncfile_path '/' files(i).name]);
    run_all_svd([ncfile_path '/' files(i).name]);
    %run_plots([ncfile_path '/' files(i).name]);
    %run_all_peaks([ncfile_path '/' files(i).name]);
end    
