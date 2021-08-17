function [ts] = ts_plot(stat,fname,lf,rf);

lfs=num2str(lf);
rfs=num2str(rf);
lfs=lfs(3:end);
rfs=rfs(3:end);

[p,~,~]=fileparts(mfilename('fullpath'));

outfile=strcat(p,'/results/',fname,'_',lfs,'-',rfs,'_ts.png');

[~,f0,~]=get_peaks(stat,lf,rf);

[N,M]=size(stat.X);

[w,c]=dpss(N,stat.p,stat.K);

%for k=2:stat.K
%    w(:,k)=(-1)*w(:,k);
%end    

%ts=real(mean(stat.Asum(f0,:),1));
ts=mean(mean(stat.ts(f0,:,:),1),3);
%ts=mean(mean(stat.signl(f0,:,:),1),3);
ts=squeeze(ts);
N=max(size(ts));
%N=70;
%eyear=2010;

fig=figure(2);
cla, clf

%plot([eyear-N:eyear],ts(end-N:end),'k','linewidth',1.5);
plot(ts,'linewidth',1.5);
ax=gca;
ax.FontSize=16;
%plot(amp(1:end));
%plot(exps(1:end));
%plot(w(:,1));
%plot(w(:,2));
%plot(w(:,3));
xlabel('Month','FontSize',14,'FontWeight','Bold');
ylabel('Amplitude','FontSize',14,'FontWeight','Bold');
xl = get(gca,'XLabel');
xlFontSize = get(xl,'FontSize');
xAX = get(gca,'XAxis');
set(xAX,'FontSize', 14)
set(xl, 'FontSize', xlFontSize);
yl = get(gca,'XLabel');
ylFontSize = get(yl,'FontSize');
yAX = get(gca,'YAxis');
set(yAX,'FontSize', 14)
set(yl, 'FontSize', ylFontSize);

title('INTERANNUAL MODE');

axis([-inf inf -inf inf]);
saveas(fig,outfile);
disp(['saved time series plot as ' outfile])

