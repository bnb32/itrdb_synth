%function MapMTMSVD(U,f0,lat,lon);
cla
%map1 = mean(U(f0-mtmstat.f0bwl:f0+mtmstat.f0bwl,:));
map1 = (U(f0,:));
absmap = abs(map1);
% map1(find(absmap <.2)) = nan;
% absmap(find(absmap <.2)) = nan;
% map1(find(absmap <.2)) = nan;

%EOFMAP
[x,y] = summermap(absmap,min(absmap),max(absmap),lat,lon);

%Plot arrows indicating phase and strength:
m_quiver(x,y,real(map1'),imag(map1'),.5,'Color','k','LineWidth',1.3)