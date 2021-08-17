clear
close all
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


%% failsafe fit curve on/off
failsafe=false;
%% run in parallel on/off
run_parallel=true;

% deterending parameters (for dtrendrw):
% 95% variance removed at 67% of sample length (conservative)
wavelength=67;
pvar=.95;

pdfit=[];

%ratio or difference detrending
%    ==1 ratio
%    ==2 difference
specs_DT_1.kopts(1)=1;
specs_DT_2.kopts(1)=1;

%empirical growth curve option
%    ==1 %N spline
%    ==8 modified negative exponential

specs_DT_1.kopts(2)=8;
specs_DT_2.kopts(2)=1;

%variance detrending
%    ==1 yes
%    ==2 no
specs_DT_1.kopts(3)=1;
specs_DT_2.kopts(3)=1;


specs_DT_1.yrpith=[];
specs_DT_2.yrpith=[];
specs_DT_1.pdfit=[];
specs_DT_2.pdfit=[];
%specs_DT_1.Splinespecs=[];
specs_DT_1.Splinespecs=[wavelength pvar];
specs_DT_2.Splinespecs=[wavelength pvar];
specs_DT_1.gcrit=1.0000e-03;
specs_DT_2.gcrit=1.0000e-03;

% minimim segment length (must be odd)
lenthresh=251;
%lenthresh=11;
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
plist=itrdb_meta.path;

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

%setup nc output filename based on detrending choice and min seg length

cpvar=num2str(pvar);
cpvar=cpvar(3:end);

if specs_DT_1.kopts(2)==1
    dtrend_type_1='spline';
elseif specs_DT_1.kopts(2)==8
    dtrend_type_1='exp';
else
    dtrend_type_1='other';
end    

if specs_DT_2.kopts(2)==1
    dtrend_type_2='spline';
elseif specs_DT_2.kopts(2)==8
    dtrend_type_2='exp';
else
    dtrend_type_2='other';
end    

dtrend_name=[dtrend_type_1 dtrend_type_2 cpvar 'v' num2str(wavelength)];
minseg_name=[num2str(lenthresh)];

ncout_name=[results_root_dir '/itrdb_Lyear_' dtrend_name '_lateHolocene_' minseg_name '_' num2str(itrdb.all_time(1)) '-' num2str(itrdb.all_time(end)) '.nc'];

% set up flags for rwl and crn failures:

import_fails=false(length(itrdb_meta.filename),1);
dtrend_fails=false(length(itrdb_meta.filename),1);
crn_fails=false(length(itrdb_meta.filename),1);

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

nr=length(flist);
%nr=10;
nl=1;

for i=nl:nr
    filename=flist{i};
    remote_path=plist{i};
    if exist(strcat(local_path,filename))~=2
        try 
            disp('------------------------------------')
            disp(['Now downloading:' filename])
            disp('------------------------------------')
	    ftp_rwl(remote_path,filename);
            pause(1.5)
	catch
            warning(['problem downloading ' filename])
	end
    end
end   

if run_parallel==true;
    delete(gcp);
    ncores=feature('numcores');
    parpool(ncores);
end

parfor i=nl:nr
%for i=1:nr
    filename=flist{i};
        disp('------------------------------------')
        disp(['Now processing: ' filename])
        disp('------------------------------------')
    try 
	rwlResults=rwl2tsm([local_path filename]);
    catch
        warning([filename ' will not work'])
        %error('STOP!')
        pause(.5)
        badrwls=[badrwls;{filename}];
        %badrwls(i)=filename;
	import_fails(i)=true;
        disp(['Import failed: ' filename])
    end

    if ~(import_fails(i));
        disp(['Import successful: ' filename])
        sname=filename;
        sname(end-3:end)=[];
        X=rwlResults.X;
	yrX=rwlResults.yrX;
        nms=rwlResults.nms;
        T=rwlResults.T;
	rwlResults_filename=[results_root_dir '/rwl/' sname '_rwl'];
	parsave(rwlResults_filename,rwlResults);
	rwl_filename{i}=rwlResults_filename;
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
            disp(['Detrending: ' filename])
	    rwiResults=dtrendrw(W,yrW,idW,yrsW,pdfit,specs_DT_1,specs_OL);
        catch
            %crnerrs=[crnerrs;{lasterr}];
	    %badchrons=[badchrons;{filename}];
	    rwierrs=[rwierrs;{lasterr}]
	    %rwierrs{i}=lasterr;
	    dtrend_fails(i)=true;
            disp(['Detrend failed: ' filename])
        end
    
        if ~(dtrend_fails(i))       
            disp(['Detrend successful: ' filename])
            rwiResults_filename=[results_root_dir '/rwi/' sname '_rwi'];
            parsave(rwiResults_filename,rwiResults)
            rwi_filename{i}=rwiResults_filename;
        
            try
                disp(['Constructing chronology: ' filename])
                crnResults=sitechron1(rwiResults.X(:,~rwiResults.explode_flag),rwiResults.yrX,minlap,crnOpt,idW(~rwiResults.explode_flag));            
            catch 
		%crnerrs{i}=lasterr;
		crnerrs=[crnerrs;{lasterr}];
		%badchrons{i}=filename;
		badchrons=[badchrons;{filename}];
		crn_fails(i)=true;                
                disp(['Chronology failed: ' filename])
            end
    
            if ~(crn_fails(i))                
                disp(['Chronology constructed: ' filename])
                crnResults_filename=[results_root_dir '/crn/' sname '_crn'];
		crn_filename{i}=crnResults_filename;
                parsave(crnResults_filename,crnResults);

                [o,qAllTime,qCronTime]=intersect(itrdb.all_time,crnResults.yry);
                if ~isempty(qAllTime)
                    
		    ntrees_mat(:,i)=parmat(itrdb.nyrs,qAllTime,crnResults.ntrees(qCronTime));
		    ncores_mat(:,i)=parmat(itrdb.nyrs,qAllTime,crnResults.ncores(qCronTime));
		    crns_raw_mat(:,i)=parmat(itrdb.nyrs,qAllTime,crnResults.y1(qCronTime));
		    crns_mat(:,i)=parmat(itrdb.nyrs,qAllTime,crnResults.y(qCronTime));
		    sss_mat(:,i)=parmat(itrdb.nyrs,qAllTime,crnResults.SSSt(qCronTime));
		    eps_mat(:,i)=parmat(itrdb.nyrs,qAllTime,crnResults.EPSt(qCronTime));
		    rbar_mat(:,i)=parmat(itrdb.nyrs,qAllTime,crnResults.rbar(qCronTime));
		end
                    
		[o,qAllTime2,qRWLTime]=intersect(itrdb.all_time,rwiResults.yrX);
                if ~isempty(qAllTime2)
                    nmssng_mat(:,i)=parmat(itrdb.nyrs,qAllTime2,rings_missing(qRWLTime));
                    fmssng_mat(:,i)=parmat(itrdb.nyrs,qAllTime2,fraction_rings_missing(qRWLTime));
		end                        
		
		sstats_rbs(i)=crnResults.rbs;
               
                if length(find(crnResults.EPSt(qCronTime)>0.85))>length(qCronTime)/3.0
		
		%if length(crnResults.y1(find(crnResults.EPSt>.85,1,'first')))==1
                    sstats_first_good_eps_yr(i)=crnResults.y1(find(crnResults.EPSt>.85,1,'first'));
                    disp(['Good fit with ' dtrend_type_1 ': ' filename])
                else
		    disp(['Bad fit with ' dtrend_type_1 ': ' filename])

                    if failsafe==true;
                        try             
                            disp(['Detrending: ' filename])
	                    rwiResults=dtrendrw(W,yrW,idW,yrsW,pdfit,specs_DT_2,specs_OL);
                        catch
                            %crnerrs=[crnerrs;{lasterr}];
	                    %badchrons=[badchrons;{filename}];
	                    rwierrs=[rwierrs;{lasterr}];
			    %rwierrs{i}=lasterr;
			    dtrend_fails(i)=true;
                            disp(['Detrend failed: ' filename])
                        end
    
                        if ~(dtrend_fails(i))       
                            disp(['Detrend successful: ' filename])
                            rwiResults_filename=[results_root_dir '/rwi/' sname '_rwi'];
                            parsave(rwiResults_filename,rwiResults)
                            rwi_filename{i}=rwiResults_filename;
                        
                            try
                                disp(['Constructing chronology: ' filename])
                                crnResults=sitechron1(rwiResults.X(:,~rwiResults.explode_flag),rwiResults.yrX,minlap,crnOpt,idW(~rwiResults.explode_flag));            
                            catch 
                                disp(['Chronology failed: ' filename])
		                crn_fails(i)=true;                
                                %crnerrs{i}=lasterr;
				crnerrs=[crnerrs;{lasterr}];
		                %badchrons{i}=filename;
				badchrons=[badchrons;{filename}];
                            end
    
                            if ~(crn_fails(i))                
                                disp(['Chronology constructed: ' filename])
                                crnResults_filename=[results_root_dir '/crn/' sname '_crn'];
	                	crn_filename{i}=crnResults_filename;
                                parsave(crnResults_filename,crnResults);

                                [o,qAllTime,qCronTime]=intersect(itrdb.all_time,crnResults.yry);
                                if ~isempty(qAllTime)
                                    
	                	    ntrees_mat(:,i)=parmat(itrdb.nyrs,qAllTime,crnResults.ntrees(qCronTime));
	                	    ncores_mat(:,i)=parmat(itrdb.nyrs,qAllTime,crnResults.ncores(qCronTime));
	                	    crns_raw_mat(:,i)=parmat(itrdb.nyrs,qAllTime,crnResults.y1(qCronTime));
	                	    crns_mat(:,i)=parmat(itrdb.nyrs,qAllTime,crnResults.y(qCronTime));
	                	    sss_mat(:,i)=parmat(itrdb.nyrs,qAllTime,crnResults.SSSt(qCronTime));
	                	    eps_mat(:,i)=parmat(itrdb.nyrs,qAllTime,crnResults.EPSt(qCronTime));
	                	    rbar_mat(:,i)=parmat(itrdb.nyrs,qAllTime,crnResults.rbar(qCronTime));
	                	end
                                    
	                	    [o,qAllTime2,qRWLTime]=intersect(itrdb.all_time,rwiResults.yrX);
                                if ~isempty(qAllTime2)
                                    nmssng_mat(:,i)=parmat(itrdb.nyrs,qAllTime2,rings_missing(qRWLTime));
                                    fmssng_mat(:,i)=parmat(itrdb.nyrs,qAllTime2,fraction_rings_missing(qRWLTime));
	                	    end                        
	                	
	                	    sstats_rbs(i)=crnResults.rbs;
                               
                                %if length(crnResults.y1(find(crnResults.EPSt>.85,1,'first')))==1
                                if length(find(crnResults.EPSt(qCronTime)>0.85))>length(qCronTime)/3.0
                                    sstats_first_good_eps_yr(i)=crnResults.y1(find(crnResults.EPSt>.85,1,'first'));
                                    disp(['Good fit with ' dtrend_type_2 ': ' filename])
                                else
	                	    disp(['Bad fit with ' dtrend_type_2 ': ' filename])
	                	end
	                    end
	                end
	            end		
		end
            end
	end
    end	
end


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
badq=false(length(itrdb_meta.filename),1);
[o,ai,bi]=(intersect(badchrons',itrdb_meta.filename));
badq(bi)=true;
itrdb.nocrn=badq;
itrdb
itrdb.results
itrdb.matrices

%disp(badchrons')
save itrdb itrdb
save badchrons badchrons
save badrwls badrwls
save crnerrs crnerrs
save rwierrs rwierrs

%output in netcdf format
ncstruct2ncfile(itrdb2ncstruct(itrdb),ncout_name);
