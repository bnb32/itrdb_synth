clf
subplot(3,1,1:2)
n = sum(~isnan(V2.tsm'));
for i = 1:size(V2.tsm,2), 
    nq(i) = find(V2.tvalid(i,1) == V2.time);
end
adj = find(nq == 294);
%subplot(3,1,1:2)
stairs(V2.time,n)
n(nq(adj(1))) = NaN;
hold on
plot(1990,0:.5:30,'k')
plot(1890,0:.5:30,'k')
plot(1850,0:.5:30,'k')
plot(1790,0:.5:30,'k')
plot(1750,0:.5:30,'k')
text(V2.tvalid(:,1),n(nq)-1,V2.id,'Fontsize',6)
text(V2.tvalid(adj,1),[24 25], char(V2.id(adj)),'Fontsize',6)
xlabel('Time (Years)')
ylabel('Number of Sites')
