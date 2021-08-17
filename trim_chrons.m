function [stats,bad_sites]=trim_chrons(itrdb);

    min_lat=30.0;
    max_lat=55.0;
    min_lon=-130.0;
    max_lon=-105.0;

    bad_sites=[];
    X_in=itrdb.chrons.data;
    lats=itrdb.lat.data;
    lons=itrdb.lon.data;
    [rows,cols]=size(X_in);
    for i=1:cols

        is_bad=gt(sum(isnan(X_in(:,i)))/rows,0.10);
	is_bad=is_bad+gt(lats(i,1),max_lat);
	is_bad=is_bad+lt(lats(i,1),min_lat);
	is_bad=is_bad+gt(lons(i,1),max_lon);
	is_bad=is_bad+lt(lons(i,1),min_lon);

	if is_bad~=0;             
	    bad_sites=[bad_sites; i];
        end
    end

    X_in(:,bad_sites)=[];
    lats(bad_sites)=[];
    lons(bad_sites)=[];
    stats.X=X_in;
    stats.lat=lats;
    stats.lon=lons;

end
