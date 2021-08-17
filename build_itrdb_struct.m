clear
close all
load itrdb_meta
load baddogs
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
    
% ouput parameters
make_rwl_plots=false;
make_rwi_plots=false;
make_crn_plots=false;

save_crns=true;
crn_dir='./';
crn_name='default';


%% STOP EDITING

ns=length(itrdb_meta.filename);
site_name_q=true(ns,1);
last_name_q=true(ns,1);
country_q=true(ns,1);
species_code_q=true(ns,1);
common_name1_q=true(ns,1);
first_year_q=true(ns,1);
last_year_q=true(ns,1);
lat_q=true(ns,1);
lon_q=true(ns,1);
elev_q=true(ns,1);               
filename_q=true(ns,1);


flist=itrdb_meta.filename;
k=1;

% set up itrdb structure for everything
itrdb.creation_date=date;
itrdb.created_by=getenv('USER');
itrdb.created_using = mfilename('fullpath') ;
meta_fnames=fieldnames(itrdb_meta);
itrdb.all_time=(-6001:2010)';
itrdb.nyrs=length(itrdb.all_time);
itrdb.nsites=length(itrdb_meta.filename);
for i =1:length(meta_fnames);
    eval(['itrdb.meta_data.' meta_fnames{i} '=itrdb_meta.' meta_fnames{i} ';'])
end
itrdb.matrices.ntrees=nan(itrdb.nyrs,itrdb.nsites);
itrdb.matrices.crns=nan(itrdb.nyrs,itrdb.nsites);
itrdb.matrices.sss=nan(itrdb.nyrs,itrdb.nsites);
itrdb.matrices.eps=nan(itrdb.nyrs,itrdb.nsites);
itrdb.matrices.nmssng=nan(itrdb.nyrs,itrdb.nsites);
itrdb.matrices.fmssng=nan(itrdb.nyrs,itrdb.nsites);


itrdb.results.rwlResults=cell(itrdb.nsites,1);
itrdb.results.rwiResults=cell(itrdb.nsites,1);
itrdb.results.crnResults=cell(itrdb.nsites,1);

% first mark "badddogs"
badq=false(length(itrdb_meta.filename),1);
[o,ai,bi]=(intersect(baddogs',itrdb_meta.filename));
badq(bi)=true;
itrdb.flags.no_rwl=badq;


badchrons={};
for i =1
    
    filename=flist{i};
    if any(strcmp(filename,baddogs))
        disp(['File: ' filename ' already flagged as problematic. Skipping...'])
    else
        [X,yrX,nms,T]=rwl2tsm([local_path filename]);        
        rwlResults.X=X;
        rwlResults.yrX=yrX;
        rwlResults.nms=nms;
        rwlResults.T=T;
        
        itrdb.results.rwlResults{i}=rwlResults;
        
        Xzero=zeros(size(X));
        Xzero(X==0)=1;
        rings_missing=sum(Xzero,2);
        fraction_rings_missing=rings_missing./sum(~isnan(X),2);
        
        q=find(sum(~isnan(X)) >= lenthresh); %
        W=X(:,q);
        yrW=yrX;
        idW=nms(q);
        yrsW=T(q,2:3);

        
        try 
            rwiResults=dtrendrw(W,yrW,idW,yrsW,pdfit,specs_DT,specs_OL);
            itrdb.results.rwiResults{i}=rwiResults;
            
            crnResults=sitechron1(rwiResults.X(:,~rwiResults.explode_flag),rwiResults.yrX,minlap,crnOpt,idW);
            itrdb.results.crnResults{i}=crnResults;
            
            [o,qyr]=intersect(itrdb.all_time,crnResults.yr);
            itrdb.matrices.ntrees(qyr,i)=crnResults.ntrees;
            itrdb.matrices.crns(qyr,i)=crnResults.y;
            itrdb.matrices.sss(qyr,i)=crnResults.SSSt;
            itrdb.matrices.eps(qyr,i)=crnResults.EPSt;
            itrdb.matrices.nmssng(qyr,i)=rings_missing;
            itrdb.matrices.fmssng(qyr,i)=fraction_rings_missing;
            
            
        catch 
            badchrons{k}=filename;
            k=k+1;
        end
    end
        
    
        
    
end

% now flag all bad files for crn building 
badq=false(length(itrdb_meta.filename),1);
[o,ai,bi]=(intersect(badchrons',itrdb_meta.filename));
badq(bi)=true;
itrdb.flags.no_crn=badq;
itrdb.flags.exlude_these=itrdb.flags.no_crn & itrdb.flags.no_rwl;
itrdb
itrdb.results
itrdb.matrices



disp(badchrons')
save itrdb itrdb
save badchrons badchrons
