function comp_all(range1,range2)

%range1='1500-2010';
%range2='1700-1850';

[p,~,~]=fileparts(mfilename('fullpath'));


files1=dir([p '/results/*' range1 '_stats.mat']);
files2=dir([p '/results/*' range2 '_stats.mat']);

for i=1:6
    file1=[p '/results/' files1(i).name];
    file2=[p '/results/' files2(i).name];
    tdiff(i)=comp_sig(file1,file2,range1,range2);
end

tdiff

end
