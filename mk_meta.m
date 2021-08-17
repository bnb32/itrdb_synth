clear
close all

%cores=feature('numcores')
%parpool(cores)

fid=fopen('./tree_pre_proc/RWLsummary_ssgupdate.txt');
k=1;
while 1
    tline=fgetl(fid);
    if ~ischar(tline)
        break,
    end
    
    q=[0 find(tline==sprintf('\t')) length(tline)+1];
    
    
    % first line is header:
    if k==1
        for i =1:length(q)-1
            fname{i}=tline(q(i)+1:q(i+1)-1);
        end
        
    else
        for i=1:length(q)-1
            eval(['raw_meta.' fname{i} '= tline(q(i)+1:q(i+1)-1);']);                            
        end
        
%        itrdb_meta.collection_id(k-1,1)=eval(raw_meta.collection_id);
%        itrdb_meta.collection_code{k-1,1}=raw_meta.collection_id;
        itrdb_meta.site_name{k-1,1}=raw_meta.site_name;
        itrdb_meta.last_name{k-1,1}=raw_meta.last_name;
        itrdb_meta.country{k-1,1}=raw_meta.country;
        itrdb_meta.species_code{k-1,1}=raw_meta.species_code;
        itrdb_meta.common_name1{k-1,1}=raw_meta.common_name1;
        itrdb_meta.first_year(k-1,1)=eval(raw_meta.first_year);
        itrdb_meta.last_year(k-1,1)=eval(raw_meta.last_year);
        
        if (raw_meta.lat_ns=='S');
            itrdb_meta.lat(k-1,1)=-(eval(raw_meta.lat_deg) + eval(raw_meta.lat_min)/60);
        elseif (raw_meta.lat_ns=='N');
            itrdb_meta.lat(k-1,1)=(eval(raw_meta.lat_deg) + eval(raw_meta.lat_min)/60);
        end
        
        if (raw_meta.lon_ew=='W')
            itrdb_meta.lon(k-1,1)=-(eval(raw_meta.lon_deg) + eval(raw_meta.lon_min)/60);
        elseif (raw_meta.lon_ew=='E')
            itrdb_meta.lon(k-1,1)=(eval(raw_meta.lon_deg) + eval(raw_meta.lon_min)/60);
        end
            
        itrdb_meta.elev_m(k-1,1)=eval(raw_meta.elev_m);
        
        itrdb_meta.path{k-1,1}=raw_meta.path;
        itrdb_meta.filename{k-1,1}=raw_meta.filename;
        
    end
    
    k=k+1;
    
end
fclose(fid)

save itrdb_meta itrdb_meta
