function [pks,locs,freqs]=get_peaks(stats,lf,rf)

    %stats=load(file);
    %[p,l]=findpeaks(stats.stats.lfv);

    pks=[];
    locs=[];
    freqs=[];
    [~,lidx]=min(abs(stats.f-lf));
    [~,ridx]=min(abs(stats.f-rf));

    %for i=1:length(l)
    %    f=stats.f(l(i));
    %    if ge(f,lf) && le(f,rf);
    %	    pks=[pks; p(i)];
    %	    locs=[locs; l(i)];
    %	    freqs=[freqs; f];
    %	end
    %end
    locs=linspace(lidx,ridx,ridx-lidx+1);
    freqs=stats.f(locs);

end    
