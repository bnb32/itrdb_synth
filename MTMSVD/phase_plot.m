function phase_plot(stats,fname,lf,rf)

[p,~,~]=fileparts(mfilename('fullpath'));

addpath([p '/../../m_map']);
addpath([p '/../']);

lfs=num2str(lf);
rfs=num2str(rf);
lfs=lfs(3:end);
rfs=rfs(3:end);

outfile=strcat(p,'/results/',fname,'_',lfs,'-',rfs,'_map.png');

[~,f0,~]=get_peaks(stats,lf,rf);

eof=stats.eof;
%[~,f0]=max(stats.lfv);

fig=figure(3);
cla,clf

map1 = sum(eof(f0,:),1);
absmap = abs(map1);

%Plot arrows indicating phase and strength:

rmap=real(map1);
imap=imag(map1);

[y,x]=meshgrid(5:-1:1);

x=reshape(x,[1,25]);
x=x(end:-1:1);
y=reshape(y,[1,25]);

quiver(x, y, rmap,imap,'Color','b','LineWidth',1.5,'AutoScaleFactor',0.75)
%m_grid('linestyle','none','tickdir','out','linewidth',1);
ax=gca;
ax.FontSize=16;
set(gca,'TickLength',[0 0]);
set(gca,'XTickLabel',[],'YTickLabel',[]);
title('INTERDECADAL MODE');

saveas(fig,outfile);
disp(['saved quiver plot as ' outfile])

end
