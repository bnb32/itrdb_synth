function [results,dt_fail,rwierr,idW,rings_missing,fraction_rings_missing,rwl_file]=dtrend_rwi(rwlResults,filename,sname,results_root_dir,imp_fail,specs_DT,specs_OL,lenthresh,pdfit);
    results.W=nan;
    results.X=nan;
    results.G=nan;
    results.yrX = nan;
    results.explode=false;
    results.explode_flag=false;
    dt_fail=true;
    rwierr='';
    rwl_file='';
    idW=nan;
    rings_missing=nan;
    fraction_rings_missing=nan;
    if ~(imp_fail);
        sname=filename;
        sname(end-3:end)=[];
        X=rwlResults.X;
	yrX=rwlResults.yrX;
        nms=rwlResults.nms;
        T=rwlResults.T;
	rwl_file=[results_root_dir '/rwl/' sname '_rwl'];
	parsave(rwl_file,rwlResults);
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
	    results=dtrendrw(W,yrW,idW,yrsW,pdfit,specs_DT,specs_OL);
	    dt_fail=false;
        catch
            %crnerrs=[crnerrs;{lasterr}];
	    %badchrons=[badchrons;{filename}];
	    rwierr=lasterr;
	    %rwierrs{i}=lasterr;
	    dt_fail=true;
        
	end    
    end
end    

