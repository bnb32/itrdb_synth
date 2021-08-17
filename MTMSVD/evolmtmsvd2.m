%Evolving MTM-SVD Spectrum. Should do the following:
%1) For any winow length, W, ends should be padded with W/2 zeros
%2) Should take steps of "step" length,
%3) Each LFV should be a measure of frequency centered at step +/- W/2

%Parameters for evolving spectrum:
W = 40;%Window size for padding, must be even
zstep =1;
padlen = round(W/2);
nfevol = N;
ufevol = round(nfevol/2);
fevol = [0:nfevol]/nfevol;
tvec = 1:zstep:N; %time vector where W is centered and MTM-SVD is performed
tvecZ = padlen +tvec; %time vector where W is centered and MTM-SVD is performed in BigZ
nest = round(N/zstep);
endopt =3; %End option: 1-Padded with zeros, 2-reflect data about endpoints, 3-truncate

%Pad ends of Data with Zeros;
if endopt ==1;   
   Padmtx = zeros(padlen, M);
elseif endopt==2;
   Padmtx = flipud(X(1:padlen,:));
elseif endopt ==3;
   clear tvec tvecZ
   Padmtx = [];
   tvec =padlen+1:zstep:N-padlen;
   tvecZ = tvec; 
   nest = round(length(tvec)/zstep);
else
    error('Invalid End Option');
end

BigZ = [Padmtx; mtmstat.Xhat; Padmtx];

%Perform iterative MTM-SVD
for i = 1:nest
    j = tvecZ(i); %Index in BigZ;
    q = j-padlen:j+padlen;
    qrec(i,:) = q;
    jrec(i,:) = j;
    h = hamming(W+1);
    [jnk, hmtx] = meshgrid(1:M,h);
    Ztemp = BigZ(q,:).*hmtx;
    lfvevol = mtmsvdf1(Ztemp,nfevol,'N');
    Fmtx(i,:) = lfvevol;
end

figure(4)
clf
subplot(3,1,1:2)
pfrayvec(1:length(tvec)) = mtmstat.p*1/W; %Vector of p x (Rayleigh frequency) for plotting
plot(time(tvec), pfrayvec,'Color',[.7 .7 .7],'LineStyle','--','LineWidth',2)
hold on
legend('f_r')
caxis([.3 .8])
pcolor(time(tvec),fevol(1:ufevol),Fmtx')
shading('interp');
xlabel('Time (years)')
ylabel('Frequency (cycles/year)')
colorbar
axis([time(tvec(1)) time(tvec(end)) 0 fevol(ufevol)])

clopt =1;
if clopt == 1
    cl = p95(end);
    contour(time(tvec),fevol(1:ufevol),Fmtx',[cl cl],'k','LineWidth',2)
else
end
plot(time(tvec), pfrayvec,'Color',[.7 .7 .7],'LineStyle','--','LineWidth',2)
title(['Evolving LFV Spectrum for ' yrttl],'FontSize',14)