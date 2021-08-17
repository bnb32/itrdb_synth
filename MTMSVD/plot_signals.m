function plot_signals(stats,outname,syear,eyear,pidx,lf,rf,dt)

    [freqs]=get_all_peaks(stats,pidx,lf,rf);

    for i=1:length(freqs)
        [R,RC,RP,~]=mtm_svd_bandrecon(stats.X,2,3,dt,freqs(i),1,0,2,0);
	if i==1
	    EOFall=RP.eof;
	    Rall=R;
	else
	    EOFall=EOFall+RP.eof;
	    Rall=Rall+R;
	end    
	
	title=strcat(['Frequency: ',num2str(freqs(i)),' (cyc/yr)']);
	ncdat.lon=stats.lon;
	ncdat.lat=stats.lat;
	fr_name=num2str(freqs(i));
	fr_name=fr_name(3:end);
	mkMTMSVDplots(RP.eof,ncdat,strcat([outname,'_phase_',fr_name,'.png']),title);
	ts_plot_gen(mean(R,1),strcat([outname,'_ts_',fr_name,'.png']),title,syear,eyear);

        %if i==length(freqs)
	%    Rall=Rall/length(freqs);
	%    EOFall=EOFall/length(freqs);
	%    mkMTMSVDplots(EOFall,ncdat,strcat([outname,'_phase_avg.png']),'Average');
	%    ts_plot_gen(mean(Rall,1),strcat([outname,'_ts_avg.png']),'Average',syear,eyear);
        %end

    end	

    

