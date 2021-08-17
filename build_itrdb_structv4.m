function build_itrdb_structv4(routine,syear,eyear)

load itrdb_meta

%% EDIT HERE
% file/path parameters:
trfunc_path='./trfuncs/';
Trunk_path='./trunk/';
no_overwrite=false;
local_path='./rwls_raw/';
results_dir_name='results';
results_root_dir=[pwd '/' results_dir_name];

addpath(trfunc_path);
addpath(Trunk_path);


if ~exist(results_root_dir)
    disp(['Creating directory: "' results_root_dir '/"'])
    mkdir(results_root_dir);

    disp(['Creating directory: "' results_root_dir '/rwl/"'])
    mkdir([results_root_dir '/rwl/']);
    
    disp(['Creating directory: "' results_root_dir '/rwi/"'])
    mkdir([results_root_dir '/rwi/']);
    
    disp(['Creating directory: "' results_root_dir '/crn/"'])
    mkdir([results_root_dir '/crn/']);
elseif exist(results_root_dir) && no_overwrite
    error(['Results directory: "' results_root_dir '/" exists. Exiting to prevent loss of data.'])
elseif exist(results_root_dir) && ~no_overwrite
     warning(['Results directory: "' results_root_dir '/" exists. Data may be over-written!!!'])   

end

%% dtrending schemes and domain variables

min_eps=0.85;

if routine==1;
    failsafe=false;
    specs_DT_1.kopts(2)=1;
    wavelength=67;
    pvar=0.5;
    lenthresh=0;
    eps_filt=false;
elseif routine==2;
    failsafe=false;
    specs_DT_1.kopts(2)=4;
    minyrs=30;
    pvar=0.95;
    lenthresh=0;
    eps_filt=false;
elseif routine==3;
    failsafe=true;
    specs_DT_1.kopts(2)=8;
    wavelength=67;
    pvar=0.5;
    lenthresh=0;
    eps_filt=false;
elseif routine==4.1;
    failsafe=false;
    specs_DT_1.kopts(2)=1;
    wavelength=100;
    pvar=0.95;
    lenthresh=251;
    eps_filt=true;
    min_eps=0.9;
elseif routine==4.2;
    failsafe=false;
    specs_DT_1.kopts(2)=4;
    minyrs=250;
    pvar=0.95;
    lenthresh=251;
    eps_filt=true;
    min_eps=0.9;
elseif routine==4.3;
    failsafe=true;
    specs_DT_1.kopts(2)=8;
    wavelength=100;
    pvar=0.95;
    lenthresh=251;
    eps_filt=true;
    min_eps=0.9;
end


%% failsafe fit curve on/off

%% run in parallel on/off
run_parallel=true;
%% verbose on/off
verbose=false;
%% filter chrons based on eps

pdfit=[];

%ratio or difference detrending
%    ==1 ratio
%    ==2 difference
specs_DT_1.kopts(1)=1;
specs_DT_2.kopts(1)=1;

%empirical growth curve option
%    ==1 %N spline
%    ==4 %nyrs spline
%    ==8 modified negative exponential

specs_DT_2.kopts(2)=1;

% deterending parameters (for dtrendrw):
% 95% variance removed at 67% of sample length (conservative)

%variance detrending
%    ==1 yes
%    ==2 no
specs_DT_1.kopts(3)=1;
specs_DT_2.kopts(3)=1;

specs_DT_1.yrpith=[];
specs_DT_2.yrpith=[];
specs_DT_1.pdfit=[];
specs_DT_2.pdfit=[];

if specs_DT_1.kopts(2)==1;
    specs_DT_1.Splinespecs=[wavelength pvar];
elseif specs_DT_1.kopts(2)==4;
    specs_DT_1.Splinespecs=[minyrs nan];
elseif specs_DT_1.kopts(2)==8;
    specs_DT_1.Splinespecs=[];
end
if failsafe && specs_DT_2.kopts(2)==1;
    specs_DT_2.Splinespecs=[wavelength pvar];
elseif failsafe && specs_DT_2.kopts(2)==4;
    specs_DT_2.Splinespecs=[minyrs nan];
end
specs_DT_1.gcrit=1.0000e-03;
specs_DT_2.gcrit=1.0000e-03;

% minimim segment length (must be odd)
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


flist=itrdb_meta.filename(1:ns);

plist=itrdb_meta.path(1:ns);

% set up itrdb structure for everything
itrdb.creation_date=date;
itrdb.created_by=getenv('USER');
itrdb.created_using = mfilename('fullpath') ;
meta_fnames=fieldnames(itrdb_meta);
itrdb.all_time=(syear:eyear)';
itrdb.nyrs=length(itrdb.all_time);
itrdb.nsites=ns;
%itrdb.nsites=length(itrdb_meta.filename);

for i =1:length(meta_fnames);
    eval(['itrdb.meta_data.' meta_fnames{i} '=itrdb_meta.' meta_fnames{i} '(1:ns);'])
end

%setup nc output filename based on detrending choice and min seg length

cpvar=num2str(pvar*100);
%cpvar=cpvar(3:end);

if specs_DT_1.kopts(2)==1 || specs_DT_1.kopts(2)==4;
    dtrend_type_1='spline';
elseif specs_DT_1.kopts(2)==8
    dtrend_type_1='exp';
else
    dtrend_type_1='other';
end    

if specs_DT_2.kopts(2)==1 || specs_DT_2.kopts(2)==4;
    dtrend_type_2='spline';
elseif specs_DT_2.kopts(2)==8
    dtrend_type_2='exp';
else
    dtrend_type_2='other';
end    

if failsafe && specs_DT_2.kopts(2)==1;
    dtrend_name=[dtrend_type_1 dtrend_type_2 cpvar 'pct' num2str(wavelength) 'wl'];
elseif failsafe && specs_DT_2.kopts(2)==4;
    dtrend_name=[dtrend_type_1 dtrend_type_2 '95pct' num2str(minyrs) 'yrs'];
elseif specs_DT_1.kopts(2)==1;
    dtrend_name=[dtrend_type_1 cpvar 'pct' num2str(wavelength) 'wl'];
elseif specs_DT_1.kopts(2)==4;
    dtrend_name=[dtrend_type_1 '95pct' num2str(minyrs) 'yrs'];
end

minseg_name=[num2str(lenthresh)];

ncout_name=[results_root_dir '/itrdb_Lyear_' dtrend_name '_lateHolocene_seglen' minseg_name '_' num2str(itrdb.all_time(1)) '-' num2str(itrdb.all_time(end)) '.nc'];

rwl_filename=cell(itrdb.nsites,1);
rwi_filename=cell(itrdb.nsites,1);
crn_filename=cell(itrdb.nsites,1);

%ntrees_mat=nan(itrdb.nyrs*itrdb.nsites,1);
ntrees_mat=nan(itrdb.nyrs,itrdb.nsites);
ncores_mat=nan(itrdb.nyrs,itrdb.nsites);
rbar_mat=nan(itrdb.nyrs,itrdb.nsites);
crns_mat=nan(itrdb.nyrs,itrdb.nsites);
crns_raw_mat=nan(itrdb.nyrs,itrdb.nsites);
eps_mat=nan(itrdb.nyrs,itrdb.nsites);
sss_mat=nan(itrdb.nyrs,itrdb.nsites);
nmssng_mat=nan(itrdb.nyrs,itrdb.nsites);
fmssng_mat=nan(itrdb.nyrs,itrdb.nsites);

sstats_rbs=nan(itrdb.nsites,1);
sstats_first_good_eps_yr=nan(itrdb.nsites,1);

%badchrons=cell(itrdb.nsites,1);
%badrwls=cell(itrdb.nsites,1);
%crnerrs=cell(itrdb.nsites,1);
%rwierrs=cell(itrdb.nsites,1);

badchrons=[];
badrwls=[];
crnerrs=[];
rwierrs=[];
rwlerrs=[];

nl=1;
nr=length(flist);

for i=nl:nr
    filename=flist{i};
    remote_path=plist{i};
    lfile=strcat(local_path,filename);
    rfile=strcat(remote_path,filename);
    if exist(lfile)~=2
        try 
	    cmd=['wget -O ' lfile ' ' rfile];
	    system(cmd);
            % ftp_rwl(remote_path,filename);
            %pause(1.5)
	    disp('------------------------------------')
            disp(['Now downloading: ' filename])
            disp('------------------------------------')
	catch
            disp('------------------------------------')
            warning(['Problem downloading: ' filename])
            disp('------------------------------------')
	end
    end
end   

if run_parallel==true && isempty(gcp('nocreate'));
    ncores=feature('numcores');
    parpool(ncores);
end

% set up flags for rwl and crn failures:

import_fails=false(length(flist),1);
dtrend_fails=false(length(flist),1);
crn_fails=false(length(flist),1);

parfor i=nl:nr
%for i=nl:nr
    filename=flist{i};
    sname=filename;
    sname(end-3:end)=[];
    disp('------------------------------------')
    disp(['Now processing: ' filename])
    disp('------------------------------------')

    [rwlResults,imp_fail,badrwl,rwlerr]=import_rwl(local_path,filename);

    import_fails(i)=imp_fail;

    rwlerrs=[rwlerrs;{rwlerr}];
    
    %{
    if imp_fail; 
        rwlerr
	return
    end
    %}

    if imp_fail==true && verbose==true;
        disp(['import failed: ' filename])
    elseif verbose==true;
        disp(['import successful: ' filename])
    end

    badrwls=[badrwls;{badrwl}];

    [rwiResults,dt_fail,rwierr,idW,rings_missing,fraction_rings_missing,rwl_file]=dtrend_rwi(rwlResults,filename,sname,results_root_dir,imp_fail,specs_DT_1,specs_OL,lenthresh,pdfit);

    dtrend_fails(i)=dt_fail;
    
    if dt_fail==true && verbose==true;;
        disp(['detrend failed with: ' filename])
    elseif verbose==true;
        disp(['detrend successful with: ' filename])
    end

    if failsafe==false;
        rwierrs=[rwierrs;{rwierr}];
        rwl_filename{i}=rwl_file;
    end

    if (import_fails(i)==false);

        [crnResults,crn_fail,crnerr,badchron,rwi_file]=build_crn(rwiResults,filename,sname,results_root_dir,dt_fail,minlap,crnOpt,idW);

        crn_fails(i)=crn_fail;
    
        if crn_fail==true && verbose==true;
            disp(['chronology failed: ' filename])
        elseif verbose==true;
            disp(['chronology successful: ' filename])
        end

        if failsafe==false;
            crnerrs=[crnerrs;{crnerr}];
            badchrons=[badchrons;{badchron}];
            rwi_filename{i}=rwi_file;
        end

        [ntrees,ncores,crns_raw,crns,sss,eps,rbar,nmssng,fmssng,rbs,first_good_eps_yr,good_fit,crn_file]=store_crn_data(crnResults,rwiResults,crn_fail,results_root_dir,sname,filename,itrdb,dtrend_type_1,rings_missing,fraction_rings_missing,eps_filt,min_eps);

        if good_fit==true && verbose==true;
            disp(['good fit with ' dtrend_type_1 ': ' filename])
        elseif verbose==true;
            disp(['bad fit with ' dtrend_type_1 ': ' filename])
        end
    end	

    if ((failsafe==true && import_fails(i)==false) && (good_fit==false || dt_fail || crn_fail));
        [rwiResults,dt_fail,rwierr,idW,rings_missing,fraction_rings_missing,rwl_file]=dtrend_rwi(rwlResults,filename,sname,results_root_dir,imp_fail,specs_DT_2,specs_OL,lenthresh,pdfit);

        dtrend_fails(i)=dt_fail;
	
	if dt_fail==true && verbose==true;
            disp(['detrend failed with: ' filename])
        elseif verbose==true;
            disp(['detrend successful with: ' filename])
        end
   
        rwierrs=[rwierrs;{rwierr}];
        rwl_filename{i}=rwl_file;

        [crnResults,crn_fail,crnerr,badchron,rwi_file]=build_crn(rwiResults,filename,sname,results_root_dir,dt_fail,minlap,crnOpt,idW);

	crn_fails(i)=crn_fail;

        if crn_fail==true && verbose==true;
            disp(['chronology failed: ' filename])
        elseif verbose==true;
            disp(['chronology successful: ' filename])
        end

        crnerrs=[crnerrs;{crnerr}];
        badchrons=[badchrons;{badchron}];
        rwi_filename{i}=rwi_file;

        [ntrees,ncores,crns_raw,crns,sss,eps,rbar,nmssng,fmssng,rbs,first_good_eps_yr,good_fit,crn_file]=store_crn_data(crnResults,rwiResults,crn_fail,results_root_dir,sname,filename,itrdb,dtrend_type_2,rings_missing,fraction_rings_missing,eps_filt,min_eps);

        if good_fit==true && verbose==true;
            disp(['good fit with ' dtrend_type_2 ': ' filename])
        elseif verbose==true;
            disp(['bad fit with ' dtrend_type_2 ': ' filename])
        end
    
    end


    if (import_fails(i)==false && dtrend_fails(i)==false && crn_fails(i)==false);

        crn_filename{i}=crn_file;
        ntrees_mat(:,i)=ntrees;
        ncores_mat(:,i)=ncores;
        crns_raw_mat(:,i)=crns_raw;
        crns_mat(:,i)=crns;
        sss_mat(:,i)=sss;
        eps_mat(:,i)=eps;
        rbar_mat(:,i)=rbar;
        nmssng_mat(:,i)=nmssng;
        fmssng_mat(:,i)=fmssng;
        sstats_rbs(i)=rbs;
        sstats_first_good_eps_yr(i)=first_good_eps_yr;
    end	
end

%if run_parallel==true && ~isempty(gcp('nocreate'));
%    delete(gcp);
%end

%fill in itrdb struct with temp variables
itrdb.flags.import_fails=import_fails;
itrdb.flags.dtrend_fails=dtrend_fails;
itrdb.flags.crn_fails=crn_fails;

itrdb.results.rwlResults_filename=rwl_filename;
itrdb.results.rwiResults_filename=rwi_filename;
itrdb.results.crnResults_filename=crn_filename;

itrdb.matrices.ntrees=ntrees_mat;
itrdb.matrices.ncores=ncores_mat;
itrdb.matrices.crns_raw=crns_raw_mat;
itrdb.matrices.crns=crns_mat;
itrdb.matrices.sss=sss_mat;
itrdb.matrices.eps=eps_mat;
itrdb.matrices.rbar=rbar_mat;
itrdb.matrices.nmssng=nmssng_mat;
itrdb.matrices.fmssng=fmssng_mat;

itrdb.site_stats.rbs=sstats_rbs;
itrdb.site_stats.first_good_eps_yr=sstats_first_good_eps_yr;

% now flag all bad files for crn building 
badq=false(length(flist),1);
[o,ai,bi]=(intersect(badchrons',itrdb_meta.filename(1:ns)));
badq(bi)=true;
itrdb.nocrn=badq;
itrdb;
itrdb.results;
itrdb.matrices;

if length(badchrons)>0
    badchrons=badchrons(find(~cellfun(@isempty,badchrons)));
end
if length(badrwls)>0
    badrwls=badrwls(find(~cellfun(@isempty,badrwls)));
end
if length(crnerrs)>0
    crnerrs=crnerrs(find(~cellfun(@isempty,crnerrs)));
end
if length(rwierrs)>0
    rwierrs=rwierrs(find(~cellfun(@isempty,rwierrs)));
end
if length(rwlerrs)>0
    rwlerrs=rwlerrs(find(~cellfun(@isempty,rwlerrs)));
end

%disp(badchrons')
%disp(['number of badchrons = ' length(badchrons)])
%itrdb.badchrons=badchrons;
%itrdb.badrwls=badrwls;
%itrdb.rwlerrs=rwlerrs;
%itrdb.crnerrs=crnerrs;
%itrdb.rwierrs=rwierrs;
save itrdb itrdb

itrdb_nc=itrdb2ncstruct(itrdb);

%output in netcdf format
ncstruct2ncfile(itrdb_nc,ncout_name);

end
