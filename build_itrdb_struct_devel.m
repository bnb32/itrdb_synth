clear
close all
load itrdb_meta
load baddogs
%% EDIT HERE
% file/path parameters:
no_overwrite=false;
local_path='./';
results_dir_name='results';
results_root_dir=[pwd '/' results_dir_name];

if ~exist(results_root_dir)
    disp(['Creating directory: "' results_root_dir '/"'])
    mkdir(results_root_dir);

    disp(['Creating directory: "' results_root_dir '/rwl/"'])
    mkdir([results_root_dir '/rwl/']);
    
    disp(['Creating directory: "' results_root_dir '/rwi/"'])
    mkdir([results_root_dir '/rwi/']);
    
    disp(['Creating directory: "' results_root_dir '/crn/"'])
    mkdir([results_root_dir '/crn/']);
end
    
if exist(results_root_dir) && no_overwrite
    error(['Results directory: "' results_root_dir '/" exists. Exiting to prevent loss of data.'])
elseif exist(results_root_dir) && ~no_overwrite
     warning(['Results directory: "' results_root_dir '/" exists. Data may be over-written!!!'])   

end


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
itrdb.matrices.nmssgn=nan(itrdb.nyrs,itrdb.nsites);
itrdb.matrices.fmssng=nan(itrdb.nyrs,itrdb.nsites);

itrdb.site_stats.rbs=nan(itrdb.nsites,1);
itrdb.site_stats.first_good_eps_yr=nan(itrdb.nsites,1);

itrdb.results.rwlResults_filename=cell(itrdb.nsites,1);
itrdb.results.rwiResults_filename=cell(itrdb.nsites,1);
itrdb.results.crnResults_filename=cell(itrdb.nsites,1);


% first mark "badddogs"
badq=false(length(itrdb_meta.filename),1);
[o,ai,bi]=(intersect(baddogs',itrdb_meta.filename));
badq(bi)=true;
itrdb.flags.rwl_fails=badq;


badchrons={};
for i =1:length(flist)
    
    filename=flist{i};
    if any(strcmp(filename,baddogs))
        disp(['File: ' filename ' already flagged as problematic. Skipping...'])
    else
        [X,yrX,nms,T]=rwl2tsm([local_path filename]);        
        sname=filename;
        sname(end-3:end)=[];
        rwlResults.X=X;
        rwlResults.yrX=yrX;
        rwlResults.nms=nms;
        rwlResults.T=T;
        
        rwlResults_filename=[results_root_dir '/rwl/' sname '_rwl'];
        save(rwlResults_filename,'rwlResults')
        
        itrdb.results.rwlResults_filename{i}=rwlResults_filename;
        
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
            
            rwiResults_filename=[results_root_dir '/rwi/' sname '_rwi'];
            save(rwiResults_filename,'rwiResults')
            itrdb.results.rwiResults_filename{i}=rwiResults_filename;
            
            crnResults=sitechron1(rwiResults.X(:,~rwiResults.explode_flag),rwiResults.yrX,minlap,crnOpt,idW);
            crnResults_filename=[results_root_dir '/crn/' sname '_crn'];            
            itrdb.results.crnResults_filename{i}=crnResults;
            save(crnResults_filename,'crnResults')
            
            [o,qyr]=intersect(itrdb.all_time,crnResults.yr);
            itrdb.matrices.ntrees(qyr,i)=crnResults.ntrees;
            itrdb.matrices.ncores(qyr,i)=crnResults.ncores;
            itrdb.matrices.crns_raw(qyr,i)=crnResults.y1;
            itrdb.matrices.crns(qyr,i)=crnResults.y;
            itrdb.matrices.sss(qyr,i)=crnResults.SSSt;
            itrdb.matrices.eps(qyr,i)=crnResults.EPSt;
            itrdb.matrices.nmssng(qyr,i)=rings_missing;
            itrdb.matrices.fmssng(qyr,i)=fraction_rings_missing;
            itrdb.matrices.rbar(qyr,i)=crnResults.rbar;
            
            itrdb.site_stats.rbs(i,1)=crnResults.rbs;
            itrdb.site_stats.first_good_eps_yr(i,1)=crnResults.y1(crnResults.EPSt>.85,1,'first');
            
            
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
itrdb.nocrn=badq;
itrdb
itrdb.results
itrdb.matrices



disp(badchrons')
save itrdb itrdb
save badchrons badchrons
