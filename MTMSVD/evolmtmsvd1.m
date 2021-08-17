%Evolving MTM-SVD Spectrum. Should do the following:
%1) For any winow length, W, ends should be padded with W/2 zeros
%2) Should take steps of "step" length,
%3) Each LFV should be a measure of frequency centered at step +/- W/2

%Parameters for evolving spectrum:
W = 40;%Window size for padding, must be even
nfevol = N;
ufevol = round(nevol/2);
zstep =5;
tvec = 1:zstep:N; %time vector where MTM-SVD is performed
padlen = round(W/2);
endopt =3; %End option: 1-Padded with zeros, 2-reflect data about endpoints, 3-truncate

%Pad ends of Data with Zeros;
if endopt ==2;   
    Padmtx = flipud(X(1:padlen,:));
else
    Padmtx = zeros(padlen, M);
end
BigZ = [Padmtx; mtmstat.Xhat; Padmtx];

%Perform iterative MTM-SVD
for k = tvec
    clear q j Z nt ns lfvevol
    q = padlen+[(k-padlen):(k+padlen)];
    l = (k-1)/zstep+1;
    recq(l,:) =q;
    reck(l,:) =k;
    recl(l,:) =l;
    
    Z = BigZ(q,:); %Matrix for analysis on this round
    [nt,ns] = size(Z);
    lfvevol = mtmsvdf1(Z,nfevol,'N');
    Fmtx(:,l) = lfvevol; %LFV spectrum centered at time(k)
end

%truncate ends if necessary:
if endopt ==1
    fgood = 1:length(tvec);
elseif endopt ==2
    fgood = 1:length(tvec);
elseif endopt ==3
    fgood = round(padlen/zstep)+1:size(Fmtx,2)-round(padlen/zstep);
else
    error('Invalid end option')
end
figure(2)
pcolor(time(tvec(fgood)), f(1:ufevol), Fmtx(:,fgood));
shading('interp');
ylabel('Frequency (Cycles/year)')
xlabel('Time (years)')
colorbar
hold on

clopt =1;
%Add confidence limits.
if clopt ==1
    [file, path] = uigetfile('*.mat','Select file with confidence limits:');
    load([path file]);
    cl90 = med_lfv(round(.1*length(med_lfv)));
    contour(time(tvec(fgood)),f(1:ufevol),Fmtx(:,fgood),[cl90 cl90],'k','LineWidth',2)
else end
