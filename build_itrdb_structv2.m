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
plist=itrdb_meta.path;
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
itrdb.matrices.ncores=nan(itrdb.nyrs,itrdb.nsites);
itrdb.matrices.rbar=nan(itrdb.nyrs,itrdb.nsites);
itrdb.matrices.crns=nan(itrdb.nyrs,itrdb.nsites);
itrdb.matrices.sss=nan(itrdb.nyrs,itrdb.nsites);
itrdb.matrices.eps=nan(itrdb.nyrs,itrdb.nsites);
itrdb.matrices.nmssng=nan(itrdb.nyrs,itrdb.nsites);
itrdb.matrices.fmssng=nan(itrdb.nyrs,itrdb.nsites);

itrdb.site_stats.rbs=nan(itrdb.nsites,1);
itrdb.site_stats.first_good_eps_yr=nan(itrdb.nsites,1);

itrdb.results.rwlResults_filename=cell(itrdb.nsites,1);
itrdb.results.rwiResults_filename=cell(itrdb.nsites,1);
itrdb.results.crnResults_filename=cell(itrdb.nsites,1);


% set up flags for rwl and crn failures:
itrdb.flags.import_fails=false(length(itrdb_meta.filename),1);
itrdb.flags.dtrend_fails=false(length(itrdb_meta.filename),1);
itrdb.flags.crn_fails=false(length(itrdb_meta.filename),1);



badchrons={};
k1=1;%rwl import
k2=1;%detrend
k3=1;%chronology
%%

%nr=length(flist);
nr=20

for i =1:nr

    filename=flist{i};
    remote_path=plist{i};
        disp('------------------------------------')
        disp(['Now processing: ' filename])
        disp('------------------------------------')
    try 
        if exist(strcat(local_path,filename))~=2
            disp('------------------------------------')
            disp(['Now downloading: ' filename])
            disp('------------------------------------')
	    ftp_rwl(remote_path,filename)
	end    
	
	%[X,yrX,nms,T]=rwl2tsm([local_path filename]);
         rwlResults=rwl2tsm([local_path filename]);
    catch
        warning('------------------------------------')
        warning([filename ' will not work'])
        warning('------------------------------------')
        %error('STOP!')
        pause(.5)
        badrwls{k1}=filename;
        k1=k1+1;
        itrdb.flags.import_fails(i)=true;
    end

    if ~(itrdb.flags.import_fails(i));
        sname=filename;
        sname(end-3:end)=[];
        X=rwlResults.X;
        yrX=rwlResults.yrX;
        nms=rwlResults.nms;
        T=rwlResults.T;
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
        catch
             crnerrs{k2}=lasterr;
             badchrons{k2}=filename;
             k2=k2+1;
             itrdb.flags.dtrend_fails(i)=true;
             
	     disp('------------------------------------')
             disp(['Detrend failed: ' filename])
             disp('------------------------------------')

        end
                      
        if ~(itrdb.flags.dtrend_fails(i))
            rwiResults_filename=[results_root_dir '/rwi/' sname '_rwi'];
            save(rwiResults_filename,'rwiResults')
            itrdb.results.rwiResults_filename{i}=rwiResults_filename;

             try
                crnResults=sitechron1(rwiResults.X(:,~rwiResults.explode_flag),rwiResults.yrX,minlap,crnOpt,idW(~rwiResults.explode_flag));            
             catch 
                  crnerrs{k3}=lasterr;
                  badchrons{k3}=filename;
                  k3=k3+1;
                  itrdb.flags.crn_fails(i)=true;
             end
            
            if ~(itrdb.flags.crn_fails(i))
                crnResults_filename=[results_root_dir '/crn/' sname '_crn'];            
                itrdb.results.crnResults_filename{i}=crnResults_filename;
                save(crnResults_filename,'crnResults')

                [o,qAllTime,qCronTime]=intersect(itrdb.all_time,crnResults.yry);
                if ~isempty(qAllTime)
                    itrdb.matrices.ntrees(qAllTime,i)=crnResults.ntrees(qCronTime);
                    itrdb.matrices.ncores(qAllTime,i)=crnResults.ncores(qCronTime);
                    itrdb.matrices.crns_raw(qAllTime,i)=crnResults.y1(qCronTime);
                    itrdb.matrices.crns(qAllTime,i)=real(crnResults.y(qCronTime));
                    itrdb.matrices.sss(qAllTime,i)=crnResults.SSSt(qCronTime);
                    itrdb.matrices.eps(qAllTime,i)=crnResults.EPSt(qCronTime);
                    itrdb.matrices.rbar(qAllTime,i)=crnResults.rbar(qCronTime);   
                    
                    [o,qAllTime2,qRWLTime]=intersect(itrdb.all_time,rwiResults.yrX);
                    itrdb.matrices.nmssng(qAllTime2,i)=rings_missing(qRWLTime);
                    itrdb.matrices.fmssng(qAllTime2,i)=fraction_rings_missing(qRWLTime);
                    
                    itrdb.site_stats.rbs(i,1)=crnResults.rbs;

                    if length(crnResults.y1(find(crnResults.EPSt>.85,1,'first')))==1
                        itrdb.site_stats.first_good_eps_yr(i,1)=crnResults.y1(find(crnResults.EPSt>.85,1,'first'));
                    end
                end


            end
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
