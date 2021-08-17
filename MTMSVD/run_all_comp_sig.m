function run_all_comp_sig

[p,~,~]=fileparts(mfilename('fullpath'));

matfile_path=[p '/results/'];

syear1=1500;
eyear1=2010;
syear2=1700;
eyear2=1850;

files1=dir([matfile_path '/*_' num2str(syear1) '-' num2str(eyear1) '_stats.mat']);
files2=dir([matfile_path '/*_' num2str(syear2) '-' num2str(eyear2) '_stats.mat']);

for i=1:length(files1)
    file1=[matfile_path '/' files1(i).name];
    file2=[matfile_path '/' files2(i).name];
    range1=[num2str(syear1) '-' num2str(eyear1)];
    range2=[num2str(syear2) '-' num2str(eyear2)];
    comp_sig(file1,file2,range1,range2);  

end    
