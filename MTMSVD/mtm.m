function [MTM] = mtm(sh,nw,k,dt,npad,sl,niter)
% function [MTM] = mtm(sh,nw,k,dt,npad,sl,niter)
%
% Function to determine MTM power spectrum
%
% Input:
%
%      sh       Time serie 
%      nw       Time-frequency bandwidth parameter [2]
%      k        Number of Slepian Tapers. Just the k=2*nw-1 are resistant
%               to spectral leakeage [3]. 
%      dt       Time interval (cycles/year)
%      npad     Set frequency range for DFT, frequency spacing, zero
%               padding, etc. npad must be a power of 2 greater than or
%               equal to the time series length (nscan). Larger values
%               provide a more interpolated frequency axis. If npad is left
%               equal to 0, npad will be 2^(nextpow2(nscan)+1).
%      sl       Significance levels: if sl is equal to 0, the default
%               levels of 0.99, 0.95, 0.90, 0.80 and 0.50 will be calculated.
%      niter    Number of permutation to perform [1000]
%
% output:
%
%      MTM      MTM power spectrum
%
%
% Example :     dt=1/12 % monthly time interval   
%               [LFV] =mtm(sh, 3, 2, dt, 0, [.99 .95 .90 .80 .50], 1000 );
%
%
% Marco Correa (marco.correa.r@gmail.com) 13.01.2016

%save A

tic
sh=sh(:);
[n] = length(sh);

% niveles de significancia
if sl==0
    sl=[.99 .95 .90 .80 .50];
end

% indice quantiles
q=fix(niter*sl)

% % Removing the mean and standard deviation
% vm=nanmean(sh); sh=sh-vm; 
% vs=nanstd(sh); sh=sh./vs; 
% clear vs vm

% rellenando faltantes con interpolacion lineal
if any(isnan(sh))
    x=1:n;
    ai=find(~isnan(sh));
    if length(ai)>1
        sh=interp1(x(ai),sh(ai), x);
    end
    sh(isnan(sh))=0;
    sh=sh(:);
end

% Slepian tapers
[psi, lam] = dpss(n,nw,k);

%[psi, lam] = dpss(n,.1,3);
%figure; plot(psi)
% padding with zeros
if npad==0; npad=2^(nextpow2(n)); end

% Some information
disp([' Slepian Tapers :        ' num2str(k) ])
disp([' Bandwidth parameter :   ' num2str(nw) ])
disp([' Padding :               ' num2str(npad) ])
disp([' Permutations :          ' num2str(niter) ])

disp(['Performing ....' ])

% dimensiones de la mitad del espectro
nf= npad/2;
% frecuencia de muestreo en npad intervalos
ddf=1./(npad*dt);
% vector de frecuencias
fr  = ([1:nf]-1')*ddf;

% Obtaining spectral stimates
psimat=[];
for a1=1:k
    psimat(:,a1)=lam(a1)*abs(fft(psi(:,a1).*sh,npad)).^2;
    %psimat(:,a1)=abs(fft(psi(:,a1),npad)).^2;
    
end
psimat=sum(psimat, 2)./sum(lam); 
psimat=mean(psimat, 2);
psimat=psimat(1:nf,:);


% 
%  figure; plot(fr, psimat(:,1)), xlim([.01 1])
% hold on; plot(fr-.0574, psimat(:,2), 'g'), xlim([.01 1])
% hold on; plot(fr-.09567, psimat(:,3), 'r'), xlim([.01 1])
% 
% S=abs(fft(sh,npad)).^2; S=S(1:nf,:);
% 
% 
% 
% 
% figure; plot(fr, psimat), xlim([.01 1])
% hold on; plot(fr, S/1000, 'k'), xlim([.01 1])

MTM.name='MultiTaper Power Spectrum'; 
MTM.timeinterval=dt;
MTM.bandwidth=nw;
MTM.tapers=k;
MTM.pading=npad;
MTM.quantil=sl;
MTM.permutations=niter;
MTM.spectrum=psimat;

% figure;
% semilogx(fr, MTM.spectrum, 'k')
% hold on; 
% semilogx(LFV_N.specdomain, LFV_N.conflevels, 'b')
% set(gca, 'xlim',[1/14 1])%,'xtick',  1./xt, 'xticklabel', xt)%, 'ylim',[0 120],'xtick',  1./xt, 'xticklabel', xt)
% 


disp([' >Local Fractional Variances : ' num2str(toc) ' s.' ])

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-----------------------------------
%    begin permutation resampling procedure
%P=ones(niter, k-1, nfmax)*nan;
%P=ones(niter, nf)*nan;
P=[];
R=(1/dt)*(1/n); %frecuencia de rayleight
dw=min(fr(end)/4, 2*nw*R); dw=dw/2; % ancho ventana

tic
for it=1:niter
    a=randperm(n); shr=sh(a);     % permutando serie original
    %figure; plot(sh, 'k'); hold on; plot(shr, 'y')
    psimat=[];
    for a1=1:k
        psimat(:,a1)=lam(a1)* abs(fft(psi(:,a1).*shr,npad)).^2;
    end
    psimat=sum(psimat, 2)./sum(lam);
    P(it,:)=psimat(1:nf)';
    
%    P=[P; psimat(1:nf, :)'];
end
disp([' >Permutations               : ' num2str(toc) ' s.' ])
%-----------------------------------

%-----------------------------------
%    now sort -- calculate 99%,95%,90%, and 50% quantiles
%    for local fractional variance as function 
%    of frequency and mode (1st mode corresponds to "local fractional
%    variance spectrum")
%
%    we need to resolve the secular band (f< p/tau where p = time-frequency
%    bandwidth product (e.g., p=2 for 2pi tapers) and tau is length of
%    data series in years since the number of spectral degrees of freedom
%    varies from k-1 to k over that band, owing to vanishing of
%    complex part of spectrum at zero frequency. Outside the band, spectral
%    degrees of freedom (complex) are k and each independent
%    frequency of the DFT
%    provides an independent sample which we can use to calculate the
%    confidence levels with a factor of roughly n/2 greater precision.
%
%    We average the estimates within and outside the secular band
%    separately

% Ordenando en cada frecuencia
P=sort(P);

% figure; pcolor(P(:, :)); shading flat
%plot(P(:,1))

% Niveles de significancia sin promediar
%E=P(q+2000,:);
E=P(q,:);

%figure; plot(fr, E', 'c')

  %vfbw=2*nw*(1/dt)*(1/n); % ancho de banda de la fecuencia espectral
R=(1/dt)*(1/n); %frecuencia de rayleight

% pomedio banda secular 
 ibsec=find(fr<2*nw*R);
%aa=mean(E(:,a),2); aa=repmat(aa,1,length(a));
bsec=E(:,ibsec);

E2=[];  E2(:,ibsec)=bsec;
a0=a(end)+1;

% pomedio banda no-secular
for i=1: length(fr)
    a=find(fr>=fr(i)-nw*R & fr<=fr(i)+nw*R);
    aa=mean(E(:,a),2); %aa=repmat(aa,1,length(a));
    E2(:,i)=aa;
end
%hold on; plot(fr, E2', 'g')

% tomando solo los mayores en la banda secular y adyacentes
a=max(bsec,E2(:,ibsec));
E2(:,ibsec)=a;

%
% % inidices banda secular y no secular
% freq_sec = nw/(n*dt);   % frecuencia secular
% ibs=find(fr<=freq_sec);  % indices banda secular
% ibns=find(fr>freq_sec);  % indices banda no-secular
% 
% % Promedio sobre la banda secular
% b=mean(P(q,ibs),2);
% E(:,ibs)= repmat(b, 1,length(ibs));
% 
% % Promedio sobre la banda no secular
% b=mean(P(q,ibns),2);
% E(:,ibns)= repmat(b, 1,length(ibns));
% %-----------------------------------

MTM.conflevels=E2;
MTM.specdomain=fr;

disp(['Done.' ])