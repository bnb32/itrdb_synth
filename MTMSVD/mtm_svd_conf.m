function [fr, evalper]=mtm_svd_conf(sh,nw,k,dt,niter,fran,sl,vw,npad)
% function [freq, evalper]=mtm_svd_conf(sh,nw,k,dt,niter,f2,sl,weight,npad)
%
% Determine the confidence levels for LFV spectrum (and for 2nd eigenvalue
% spectrum also based on bootstrap resampling calculations.
%
% Input:
%
%      sh       Time series matrix [rows=time columns=space]
%      nw       Time-frequency bandwidth parameter [2]
%      k        Number of Slepian Tapers. Just the k=2*nw-1 are resistant
%               to spectral leakeage [3].
%      dt       Time interval (cycles/year). Instead of the time interval
%               parameter, a regular date vector in Matlab format (as the
%               produced by datenum.m) can be entered.
%      niter    Number of permutation to perform [1000]
%      fran     fran=[f1 f2], range of frequencies (in cyc/year) to
%               evaluate. f1/f2 are the lowest/higest frequencies. If fran
%               is set to 0, the range will be set from the secular to the
%               nyquist sampling frequency: fran=[0 0.5/dt].  
%      sl       Significance levels: if sl is equal to 0, the default
%               levels of 0.99, 0.95, 0.90, 0.80 and 0.50 will be calculated.
%      vw       Weights on dataseries: the user may decide to set e.g., 
%               areal(norm=cos(latitude) [space]) weights on the constituent
%               dataseries. Otherwise uniform weights are set if weight=0
%      npad     Set frequency range for DFT, frequency spacing, zero
%               padding, etc. npad must be a power of 2 greater than or
%               equal to the time series length (nscan). Larger values
%               provide a more interpolated frequency axis. If npad is left
%               equal to 0, npad will be 2^(nextpow2(nscan)+1).
%
% output:
%   evalper     Matrix containing significance levels for first eigenvalue
%               LFV spectrum Columns contain the quantiles [sl] for the
%               local fractional variance spectrum as function of the
%               frequency.
%   fr          Frequency domain vector
%
%-------------------------------------------------------------------------- 
% A description of the theoretical basis of the MTM-SVD toolbox and some
% implementation details can be found in:
%
%  Correa-Ramirez, M. & S. Hormazabal, 2012. "MultiTaper Method-Singular
%  Value Decomposition (MTM-SVD): Variabilidad espacio–frecuencia de las
%  fluctuaciones del nivel del mar en el Pacífico Suroriental". Lat. Am. J.
%  Aquat. Res., 40[4]: 1039-1060 . 
%
% Citation of this article would be appreciated if you use the toolbox
%
% A demonstration program is included which tests many of the capabilities
% of this toolbox. To run the demonstration, type 'test_mtm t_svd'.
%
% Copyright (C) 2012, Marco Correa-Ramirez and Samuel Hormazabal
% Pontificia Universidad Catolica de Valparaiso
% Escuela de Ciencias del Mar, Valparaiso, Chile
%
% This software may be used, copied, or redistributed as long as it is not
% sold and this copyright notice is reproduced on each copy made. This
% routine is provided as is without any express or implied warranties
% whatsoever.
%
% ------------------------------------------------------------------------
% The functions of this toolbox are based on the MTM-SVD FORTRAN functions
% developed by Michael Mann and Jeffrey Park
% (http://holocene.meteo.psu.edu/Mann/tools/MTM-SVD/)
% ------------------------------------------------------------------------
%
%
% Questions or comments to:
% M. Correa-Ramirez (marco.correa.r@gmail.com)
% 
% Last update (version)
% 2014.04.02  (ver. 2.0)

save A

% dimensiones de la matiz        
[n, p]=size(sh);

% intervalo de tiempo
if length(dt)>1
    lt=(dt(end)-dt(1))/365.25;  % longitud de la serie en años
    dt=lt/length(dt);  
end

% Remueve la media y divide por la desviacion estandar
vm=nanmean(sh); vm=repmat(vm,n,1);
sh=sh-vm; clear vm
vs=nanstd(sh); vs=repmat(vs,n,1);
sh=sh./vs; clear vs

% niveles de significancia
if sl==0
    sl=[.99 .95 .90 .80 .50];
end

% indice de los cuantiles
q=fix(niter*sl);

% inicializa las series de tiempo con pesos uniformes. Se puede usar otra
% ponderacion eg. ponderacion por area [cos(latitude)] sobre cada serie
% de tiempo.

 if vw==0; vw=ones(1,p); end

% ponderando por los pesos
  W = repmat(vw(:)',n,1);
  sh = sh.*W; clear W
%-----------------------------------------------------
% Definicion de parametros
 f1 = 0;                  % minima frecuencia para el analisis
% if f2==0; f2=0.5/dt ;end % maxima frecuencia para el analisis

% f2 puede variar hasta al frecuencia de nyquist (0.5/dt). Esto permite
% ahorrar mucho tiempo de computo

 if npad==0; npad=2^(nextpow2(n)); end % padding con zeros
 ddf = 1./(npad*dt); % frecuencia de muestreo en npad intervalos
 nf = npad/2; % dimensiones de la mitad del espectro
 fr = ([1:nf]-1')*ddf;% frequecies
 
 if length(fran)==1 & fran==0
     fran=[0 0.5/dt];
 end
 if length(fran)==1 & fran~=0
     fran=[0 fran];
 end

 nfmin = near(fr,fran(1));  
 nfmax = near(fr,fran(2));% indice de corte   
 fr = fr(nfmin:nfmax);% frecuencia
%ssdfggf
%------------------------------------------------- 
% rellena datos faltantes con interpolacion lineal
if any(isnan(sh(:)))
    x=1:n;
    for i=1:p;
        t=sh(:,i);
        ai=find(~isnan(t));
        if ~isempty(ai)
        sh(:,i)=interp1(x(ai),t(ai), x)';
        end
    end
    sh(isnan(sh))=0;
end

%----------------------------------------------
% construye tapers apropiados para el intervalo de tiempo y ancho de banda
[psi] = dpss(n,nw, k);

% Redimensionando tapers
for i1=1:3
    psi2(:,:,i1)=repmat(psi(:,i1),1,p);
end
%----------------------------------------------
%  Realiza el remuestreo y permutacion

 partvar=ones(niter, nfmax-nfmin+1)*nan;
 rng('shuffle') 
for it=1:niter    
    tic
    a=randperm(n); shr=sh(a,:); % desordena en el tiempo la serie original
    
% Obtiene las estimaciones espectrales para cada uno de los tapers
% ("k" Slepian tapers)
    nev=[];
    for i1=1:k
        %psimat=repmat(psi(:,i1),1,p).*shr;   % Redimensionando tapers
        psimat=psi2(:,:,i1).*shr;        % Aplicando tapers
        psimat=fft(psimat,npad);        % Estimativos spectrales con tapers
        psimat(nfmax+1:end, :)=[];      % Conservando mitad del espectro
        skn=['sk' num2str(i1)];         % Nombre matriz de almacenamiento
        eval([skn '= psimat;']);        % Matriz
        nev=[nev skn '(j1, :)''  '];    % matriz de evalaucion
    end
    clear psimat shr skn

    % perform SVD for each frequency
    %nfmin:nfmax
    a=0;
    for j1=nfmin:nfmax
        eval(['M=[' nev '];'])
        [S] = svd(M, 0);
        %asdfgf
        a=a+1;
        partvar(it, a)=(S(1).^2)/sum(S(1:end).^2);
    end
    %sdffg
    clear sk*
    clc
    disp(['CONFIDENCE LEVELS'])
    disp(['MTM-SVD LFV SPECTRUM'])
    disp([' Slepian Tapers :        ' num2str(k) ])
    disp([' Bandwidth parameter :   ' num2str(nw) ])
    disp([' Padding :               ' num2str(npad) ])
    disp([' Total Permutations :    ' num2str(niter) ])
    disp(['Performing ....' ])
    disp(['>Permutation #: ' num2str(it) ' | Time to complete: ' sprintf('%4.1f', (niter-it)*toc/60) 'min'])
end

%-----------------------------------
disp(['Done.' ])
% load B partvar fr n nw dt q
partvar=sort(partvar);              % Ordenando en cada frecuencia
% banda secular
freq_sec = nw/(n*dt);               % frecuencia secular
ibs=near(fr, freq_sec); ibs=1:ibs;  % indices banda secular
ibns=ibs(end)+1:length(fr);         % indices banda no-secular
fray = 1/(n*dt);                    % frecuencia de Rayleigh
fbw=2*nw*fray;                      % ancho de banda de la estimacion espectral
ifbw=round(fbw/diff(fr([1 2])));    % No. indices ancho de banda

warning off
evalper=[];
for i=1:length(q)
    y=partvar(q(i),:);
    y1=medfilt1(y,ifbw);            % filtrando por hacho de banda
    % hold on; plot(fr,y1, 'r')
    y2(ibs)=mean(y(ibs));           % Promedio sobre la banda secular
    [a,d]=polyfit(fr(ibns),y1(ibns),10);% Fijando polinomio
    y2(ibns)=polyval(a,fr(ibns));
    %hold on; plot(fr,y2, 'g','LineWidth',2)
    evalper(i,:)=y2;
end
warning on


%figure; plot(fr, evalper,'LineWidth',2)

