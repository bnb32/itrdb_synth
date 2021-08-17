function phase_plot_gen(map1,outfile,title_name)

[p,~,~]=fileparts(mfilename('fullpath'));

addpath([p '/../../m_map']);
addpath([p '/../']);

fig=figure(3);
cla,clf

absmap = abs(map1);

%Plot arrows indicating phase and strength:

rmap=real(map1);
imap=imag(map1);

%[y,x]=meshgrid(1:5);
[y,x]=meshgrid(5:-1:1);

x=reshape(x,[25,1]);
%x=x(end:-1:1);
y=reshape(y,[25,1]);

quiver(x, y, rmap,imap,'Color','b','LineWidth',1.5,'AutoScaleFactor',0.75)
%m_grid('linestyle','none','tickdir','out','linewidth',1);
ax=gca;
ax.FontSize=16;
set(gca,'TickLength',[0 0]);
set(gca,'XTickLabel',[],'YTickLabel',[]);
axis([0.5 6 0.5 6]);
title(title_name);

saveas(fig,outfile);
disp(['saved quiver plot as ' outfile])

end
