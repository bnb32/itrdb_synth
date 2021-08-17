clear
close all
load badchrons

%% EDIT HERE
% file/path parameters:
local_path='./';

% space/time/network selection parameters:
tdomain=[-inf inf];
lonE=-110;
lonw=-115;
latS=35;
latN=40;

% deterending parameters (for dtrendrw):
% 95% variance removed at 67% of sample length (conservative)
wavelength=67;
pvar=.95;

pdfit=[];

specs_DT.kopts=[1; 1; 1];
specs_DT.yrpith=[];
specs_DT.pdfit=[];
specs_DT.Splinespecs=[wavelength pvar];
specs_DT.gcrit=1.0000e-03;

% minimim segment length (must be odd)
lenthresh=11;
specs_OL.myr=lenthresh-2;
specs_OL.nstdev=[3 3];
specs_OL.kremove=0;

% Chronology-building parameters:
minlap=3; %minimum overlap

% variance stabilization for time varying sample size
%        ==1 Yes 
%        ==2 No
crnOpt(1)=1;


% variance equalization of core indices before making site chronology
%         ==1 Yes
%         ==2 No
crnOpt(2)=1;
    

% QC parameters

% ouput parameters
make_rwl_plots=false;
make_rwi_plots=false;
make_crn_plots=false;

save_crns=true;
crn_dir='./';
crn_name='default';





%% STOP EDITING



flist=badchrons;
k=1;

for i =1:length(flist)
    
    filename=flist{i};
    [X,yrX,nms,T]=rwl2tsm([local_path filename]);        
    q=find(sum(~isnan(X)) >= lenthresh); %
    W=X(:,q);
    yrW=yrX;
    idW=nms(q);
    yrsW=T(q,2:3);
    
    try 
        rwlResults=dtrendrw(W,yrW,idW,yrsW,pdfit,specs_DT,specs_OL); 
        crnResults=sitechron1(rwlResults.X,rwlResults.yrX,minlap,crnOpt,idW);
    catch 
        badchrons{k}=filename;
        k=k+1;
    end
        
        
    
end

disp(badchrons)