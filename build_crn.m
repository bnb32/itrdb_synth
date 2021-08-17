function [results,crn_fail,crnerr,badchron,rwi_file]=build_crn(rwiResults,filename,sname,results_root_dir,dt_fail,minlap,crnOpt,idW);
   results.y=nan;
   results.yry=nan;
   results.y1=nan;
   results.yry1 = nan;
   results.ylo=nan;
   results.yhi=nan;
   results.ys=nan;
   results.SSSt = nan;
   results.EPSt = nan;
   results.rbar=nan;
   results.ntrees=nan;
   results.ncores=nan;
   results.rbs = nan;
   results.SSS=nan;
   results.EPS = nan;
   results.f = nan;
   crn_fail=true;
   crnerr='';
   rwi_file='';
   badchron='';
   if ~(dt_fail)       
        rwi_file=[results_root_dir '/rwi/' sname '_rwi'];
        parsave(rwi_file,rwiResults)
    
        try
            results=sitechron1(rwiResults.X(:,~rwiResults.explode_flag),rwiResults.yrX,minlap,crnOpt,idW(~rwiResults.explode_flag));            
            crn_fail=false;
	catch 
    	%crnerrs{i}=lasterr;
    	    crnerr=lasterr;
    	%badchrons{i}=filename;
    	    badchron=filename;
    	    crn_fail=true;                
        end
    end
end    
