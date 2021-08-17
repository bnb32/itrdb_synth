function download_itrdb_data

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

for i=1:length(flist)
    filename=flist{i};
    remote_path=plist{i};
    lfile=strcat(local_path,filename);
    rfile=strcat(remote_path,filename);
    if exist(lfile)~=2
        try 
            %ftp_rwl(remote_path,filename);
	    cmd=['wget -O ' lfile ' ' rfile];
	    system(cmd);
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

