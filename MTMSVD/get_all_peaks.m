function [freqs,pvals]=get_all_peaks(stats,pidx,lf,rf)
    pks=[];
    fqs=[];
    locs=[];
    pvs=[];
    [p,l]=findpeaks(stats.lfv);
    f=stats.f;
    lfv=stats.lfv;
    clv=stats.clv;
    for i=1:length(p)
        if ge(lfv(l(i)),clv(l(i),pidx)) && ge(f(l(i)),lf) && le(f(l(i)),rf);
	    fqs=[fqs;f(l(i))];
            if ge(lfv(l(i)),clv(l(i),1)); pvs=[pvs;0.99];
	    elseif ge(lfv(l(i)),clv(l(i),2)); pvs=[pvs;0.95];
	    elseif ge(lfv(l(i)),clv(l(i),3)); pvs=[pvs;0.90];
	    elseif ge(lfv(l(i)),clv(l(i),4)); pvs=[pvs;0.80];
	    elseif ge(lfv(l(i)),clv(l(i),5)); pvs=[pvs;0.50];
	    else pvs=[pvs;0.0];
	    end
	    
	end
    end	
    pvals=pvs;
    freqs=fqs;
end    
