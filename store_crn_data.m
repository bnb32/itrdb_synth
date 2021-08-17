function [ntrees,ncores,crns_raw,crns,sss,eps,rbar,nmssng,fmssng,rbs,first_good_eps_yr,good_fit,crn_file]=store_crn_data(crnResults,rwiResults,crn_fail,results_root_dir,sname,filename,itrdb,dtrend_type,rings_missing,fraction_rings_missing,eps_filt,min_eps);
    nyrs=itrdb.nyrs;
    all_time=itrdb.all_time;
    ntrees=nan(nyrs,1);
    ncores=nan(nyrs,1);
    crns_raw=nan(nyrs,1);
    crns=nan(nyrs,1);
    sss=nan(nyrs,1);
    eps=nan(nyrs,1);
    rbar=nan(nyrs,1);
    nmssng=nan(nyrs,1);
    fmssng=nan(nyrs,1);
    rbs=nan;
    first_good_eps_yr=nan;
    good_fit=false;
    crn_file='';
    if ~(crn_fail)                
        crn_file=[results_root_dir '/crn/' sname '_crn'];
        parsave(crn_file,crnResults);
             
        [o,qAllTime,qCronTime]=intersect(all_time,crnResults.yry);
        indices=find(crnResults.EPSt(qCronTime)>min_eps);
	if length(indices)>length(qCronTime)/3.0
            if ~isempty(qAllTime)
	   
                ntrees=parmat(nyrs,qAllTime,crnResults.ntrees(qCronTime));
                ncores=parmat(nyrs,qAllTime,crnResults.ncores(qCronTime));
                crns_raw=parmat(nyrs,qAllTime,crnResults.y1(qCronTime));
		sss=parmat(nyrs,qAllTime,crnResults.SSSt(qCronTime));
                eps=parmat(nyrs,qAllTime,crnResults.EPSt(qCronTime));
                rbar=parmat(nyrs,qAllTime,crnResults.rbar(qCronTime));    
	        
		if ~(eps_filt)
                    crns=parmat(nyrs,qAllTime,crnResults.y(qCronTime));
            
	        elseif eps_filt
                    tmp=crnResults.y(qCronTime);
		    crns=parmat(nyrs,qAllTime(indices),tmp(indices));
	        end 

            end
            
            [o,qAllTime2,qRWLTime]=intersect(all_time,rwiResults.yrX);
            if ~isempty(qAllTime2)
                nmssng=parmat(nyrs,qAllTime2,rings_missing(qRWLTime));
                fmssng=parmat(nyrs,qAllTime2,fraction_rings_missing(qRWLTime));
            end                        
        
            rbs=crnResults.rbs;

            %length(crnResults.y1(find(crnResults.EPSt>.85,1,'first')))==1
	    good_yrs=crnResults.y1(indices);
            sstats_first_good_eps_yr=good_yrs(1);
	    good_fit=true;
        end
    end	
end
