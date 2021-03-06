function [hidx, ctree] = mtm_svd_equvar(sh,nw,k,dt,fqr,npad,cutoff)
% function [idvar, ctree] = mtm_svd_equvar(sh,k,nw,dt,weight,npad);
%
% Determine equal variance zones from a set of multi taper frequency
% espectra distributed in a large region.  
%
% Input:
%
%      sh       Time series matrix [rows=time columns=space]
%      nw       Time-frequency bandwidth parameter [2]
%      k        Number of Slepian Tapers. Just the k=2*nw-1 are resistant
%               to spectral leakeage [3]. 
%      dt       Time interval (cycles/year)
%      fqr      Vector of frequency range: [fq(lower) fq(higher)]. If input
%               fqr=0, fqr is set from 0 to nyquist sampling frequency:
%               fqr=[0 0.5/dt].
%      npad     Set frequency range for DFT, frequency spacing, zero
%               padding, etc. npad must be a power of 2 greater than or
%               equal to the time series length (nscan). Larger values
%               provide a more interpolated frequency axis. If npad is left
%               equal to 0, npad will set to 2^(nextpow2(nscan)+1).
%      cutoff   Threshold distance for cutting. Distance is measured like
%               1-correlation coeficient, so the cutoff treshold varies
%               betwen 0-1. If cutoff=0, cutoff is set to 0.5.  
%
% 
% output:
%
%      idvar    Vector containig the clustering indices
%      ctree    Matrix generated by the clustering which is useful to
%               generates a dendrogram plot. 
%
% Marco Correa (mcorrea@dgeo.udec.cl) 01.12.2011

disp(['MTM-SVD EQUAL VARIANCE ZONES'])
disp(['Preprocesing ...         ' ])
[n,p] = size(sh);

% Removing the mean and standard deviation
vm=nanmean(sh); vm=repmat(vm,n,1);
sh=sh-vm; clear vm
vs=nanstd(sh); vs=repmat(vs,n,1);
sh=sh./vs; clear vs

% %    initialize uniform weights on dataseries -- the user
% %    may decide to set e.g., areal (ie, cos(latitude) ) weights
% %    on the constituent dataseries, etc.
% if vw==0; vw=ones(1,p); end
% 
% % ponderando por los pesos
% W=repmat(vw(:)',n,1);
% sh=sh.*W; clear W

% rellenando faltantes con interpolacion lineal
if any(isnan(sh(:)))
    x=1:n;
    for i=1:p;
        t=sh(:,i);
        ai=find(~isnan(t));
        if length(ai)>1
            sh(:,i)=interp1(x(ai),t(ai), x)';
        end
    end
    sh(isnan(sh))=0;
end

% Slepian tapers
[psi] = dpss(n,nw,k);

% padding with zeros
if npad==0; npad=2^(nextpow2(n)); end

% Some information
disp([' Slepian Tapers :        ' num2str(k) ])
disp([' Bandwidth parameter :   ' num2str(nw) ])
disp([' Padding :               ' num2str(npad) ])
disp(['Performing ....' ])

% dimensiones de la mitad del espectro
nf= npad/2;
% frecuencia de muestreo en npad intervalos
ddf=1./(npad*dt);
% vector de frecuencias
fr  = ([1:nf]-1')*ddf;

% Obtaining spectral stimates
Y=[];
for i1=1:k
    tic
    psimat=repmat(psi(:,i1),1,p);   % Redimensionando tapers
    psimat=fft(psimat.*sh,npad);    % Estimativos spectrales con tapers
    psimat(nf+1:end, :)=[];         % Conservando mitad del espectro
    psimat=abs(psimat).^2;
    Y=sum(cat(3, Y, psimat), 3);
    disp([' >Spectral stimated ' num2str(i1) ': ' num2str(toc) ' s.'])
end
clear psimat % saving memory
Y=Y./k;
disp(['Computing distances ....' ]);
tic
D = pdist(Y','correlation');
disp([' >Performed ' num2str(p*(p-1)/2) ' distances in ' num2str(toc) ' s.'])
disp(['Computing the agglomerative hierarchical structure ....' ]);
ctree = linkage(D,'complete');
% indice de grupos
hidx = cluster(ctree,'criterion','distance','cutoff',cutoff);

disp([' >Completed.'])
disp(['Clusters information:' ]);
disp([' Threshold for cutting :  ' num2str(cutoff) ])
disp([' No. of zones formed :    ' num2str(max(hidx)) ])
disp(['Done.' ])


