function [stats]=LME_prep(data)

lat=data.lat.data;
lon=data.lon.data;
time=data.time.data;
vals=data.H2OSOI.data;

min_lon=-130;
max_lon=-100;
min_lat=25.0;
max_lat=60.0;

%min_lon=-180;
%max_lon=180;
%min_lat=-90;
%max_lat=90;

for i=1:length(lon)
    if ge(lon(i),180.0)
        lon(i)=lon(i)-360.0;
    end
end

count=1;
for i=1:length(lon)
    for j=1:length(lat)
        if ge(lon(i),min_lon) && le(lon(i),max_lon) && ge(lat(j),min_lat) && le(lat(j),max_lat)
	    Xdat(:,count)=sum(vals(i,j,:,:),3);
	    lats(count)=lat(j);
	    lons(count)=lon(i);
	    count=count+1;
	end
    end
end    

stats.lat=lats;
stats.lon=lons;
stats.X=Xdat;
