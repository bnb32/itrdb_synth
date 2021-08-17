function ncstruct=itrdb2ncstruct(itrdb)

    ncstruct.dims.nsites=itrdb.nsites;
    ncstruct.dims.nyears=itrdb.nyrs;
    
    %{
    ncstruct.country.data=cellstr(itrdb.meta_data.country);
    ncstruct.country.dim_names=cellstr('nsites');
    
    ncstruct.site_name.data=cellstr(itrdb.meta_data.site_name);
    ncstruct.site_name.dim_names=cellstr('nsites');
    
    ncstruct.tree_name.data=cellstr(itrdb.meta_data.common_name1);
    ncstruct.tree_name.dim_names=cellstr('nsites');
    
    ncstruct.species_code.data=cellstr(itrdb.meta_data.species_code);
    ncstruct.species_code.dim_names=cellstr('nsites');
    
    ncstruct.raw_file_name.data=cellstr(itrdb.meta_data.filename);
    ncstruct.raw_file_name.dim_names=cellstr('nsites');
    
    ncstruct.raw_file_path.data=cellstr(itrdb.meta_data.path);
    ncstruct.raw_file_path.dim_names=cellstr('nsites');
    %}
    
    ncstruct.time.data=itrdb.all_time;
    ncstruct.time.dim_names=cellstr('nyears');
    
    ncstruct.nochron.data=uint8(itrdb.nocrn);
    ncstruct.nochron.dim_names=cellstr('nsites');
    
    ncstruct.lat.data=itrdb.meta_data.lat;
    ncstruct.lat.dim_names=cellstr('nsites');
    
    ncstruct.lon.data=itrdb.meta_data.lon;
    ncstruct.lon.dim_names=cellstr('nsites');
    
    ncstruct.elev.data=itrdb.meta_data.elev_m;
    ncstruct.elev.dim_names=cellstr('nsites');
    
    ncstruct.first_year.data=itrdb.meta_data.first_year;
    ncstruct.first_year.dim_names=cellstr('nsites');
    
    ncstruct.last_year.data=itrdb.meta_data.last_year;
    ncstruct.last_year.dim_names=cellstr('nsites');
    
    
    ncstruct.detrend_fail.data=uint8(itrdb.flags.dtrend_fails);
    ncstruct.detrend_fail.dim_names=cellstr('nsites');
    
    ncstruct.chron_fail.data=uint8(itrdb.flags.crn_fails);
    ncstruct.chron_fail.dim_names=cellstr('nsites');
    
    ncstruct.chrons.data=itrdb.matrices.crns;
    ncstruct.chrons.dim_names={'nyears','nsites'};
    
    ncstruct.chrons_raw.data=itrdb.matrices.crns_raw;
    ncstruct.chrons_raw.dim_names={'nyears','nsites'};
    
    ncstruct.eps.data=itrdb.matrices.eps;
    ncstruct.eps.dim_names={'nyears','nsites'};
    
    ncstruct.fmssng.data=itrdb.matrices.fmssng;
    ncstruct.fmssng.dim_names={'nyears','nsites'};
    
    ncstruct.ncores.data=itrdb.matrices.ncores;
    ncstruct.ncores.dim_names={'nyears','nsites'};
    
    ncstruct.nmssng.data=itrdb.matrices.nmssng;
    ncstruct.nmssng.dim_names={'nyears','nsites'};
    
    ncstruct.ntrees.data=itrdb.matrices.ntrees;
    ncstruct.ntrees.dim_names={'nyears','nsites'};
    
    ncstruct.rbar.data=itrdb.matrices.rbar;
    ncstruct.rbar.dim_names={'nyears','nsites'};
    
    ncstruct.sss.data=itrdb.matrices.sss;
    ncstruct.sss.dim_names={'nyears','nsites'};

    %ncstruct.badchrons.data=itrdb.badchrons;
    %ncstruct.badrwls.data=itrdb.badrwls;
    %ncstruct.rwlerrs.data=itrdb.rwlerrs;
    %ncstruct.crnerrs.data=itrdb.crnerrs;
    %ncstruct.rwierrs.data=itrdb.rwierrs;
    
    ncstruct.global_attributes.history=['File created by Brandon Benton']

end
