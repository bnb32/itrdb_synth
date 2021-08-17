%load ../itrdb
%load stats

function mkMTMSVDplots(map1,ncdat,outname,plot_title)

[p,~,~]=fileparts(mfilename('fullpath'));

addpath([p '/../../m_map']);
addpath([p '/../']);

lat = ncdat.lat;
lon = ncdat.lon;
nlat = length(lat);
nlon = length(lon);

lonw=-130;
lone=-100;
lats=30.0;
latn=55.0;

fig=figure(1);
cla,clf

absmap = abs(map1);
arrow_mean=mean(absmap);
%map1(find(absmap<arrow_mean))=0;

%Plot arrows indicating phase and strength:

m_proj('Mercator','longitude',[lonw lone],'latitude',[lats latn]);  
m_coast('patch', [.9 .9 .9], 'edgecolor', 'k','LineWidth',2);
hold on

sv=1;
rmap=real(map1');
imap=imag(map1');
%m_quiver(longlon,longlat,real(map1'),imag(map1'),1.5,'Color','k','LineWidth',1.3)
m_quiver(lon(1:sv:end),lat(1:sv:end),rmap(1:sv:end),imap(1:sv:end),3.5,'Color','b','LineWidth',1.0,'AutoScaleFactor',0.5)


ax=gca;
ax.FontSize=16;
title(plot_title);

%Z=diag(absmap);
%m_contour(lon,lat,Z);
m_grid('linestyle','none','tickdir','out','linewidth',1);
saveas(fig,outname);
disp(['saved quiver plot as ' outname])
end
