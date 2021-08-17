function [LFV] = mtm_svd_lfv(sh,nw,k,dt,vw,npad)
% function [LFV] = mtm_svd_lfv2(sh,k,nw,dt,weight,npad);
%
% Determine the local fractional variance spectrum LFV 
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
%      weight   Weights on dataseries: the user may decide to set e.g.,
%               areal(norm=cos(latitude) [space]) weights on the constituent
%               dataseries.  Otherwise uniform weights are set if
%               weight=0
%      npad     Set frequency range for the spectrum domain, frequency
%               spacing, zero padding, etc. npad must be a power of 2
%               greater than or equal to the time series length (nscan).
%               Larger values provide a more interpolated frequency axis.
%               If npad is left equal to 0, npad will be 2^(nextpow2(nscan)+1).
%
% output:
%
%   LFV                     Structure array containig:
%   LFV.timeinterval        Time interval (cycles/year)
%   LFV.bandwidth           Bandwidth parameter used
%   LFV.tapers              Slepian Tapers used
%   LFV.pading              Pading in the spectrum domain
%   LFV.varfreqbandwidth    Frequency variable bandwidth. Effective spectral
%                           resolution (= 2 x nw x (1/dt) x (1/n)) arround
%                           each frequency in the spectral domain.  
%   LFV.spectrum            Local fractional variance spectrum [1 x npad]
%   LFV.specdomain          frequency domain vector [1 x npad]
% 
% Example:
%
%       [LFV] = mtm_svd_lfv2(sh,2,3,.09,0,512);
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
% of this toolbox. To run the demonstration, type 'test_mtmt_svd'.
%
% Copyright (C) 2012, Marco Correa-Ramirez and Samuel Hormazabal
% Pontificia Universidad Catolica de Valparaiso
% Escuela de Ciencias del Mar, Valparaiso, Chile
%
% This software may be used, copied, or redistributed as long as it is not
% sold and this copyright notice is reproducebandd on each copy made. This
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
% 2014.04.02  (ver. 1.0)

% save A
% 
% clear all; load A

disp(['MTM-SVD LOCAL FRACTIONAL VARIANCE SPECTRUM'])
disp(['Preprocesing ...         ' ])
[n,p] = size(sh);

% Intervalo de tiempo
if length(dt)>1
    lt=(dt(end)-dt(1))/365.25;    % longitud de la serie en años
    dt=lt/length(dt);  
end

% Remueve la media y divide por la desviacion estandar
 vm = nanmean(sh); vm=repmat(vm,n,1);
 sh = sh-vm; clear vm
 vs = nanstd(sh); vs=repmat(vs,n,1);
 sh = sh./vs; clear vs

% inicializa las series de tiempo con pesos uniformes. Se puede usar otra
% ponderacion eg. ponderacion por area [cos(latitude)] sobre cada serie
% de tiempo.
  if vw==0; vw=ones(1,p); end
  
% ponderando por los pesos
  W = repmat(vw(:)',n,1);
  whos W
  sh = sh.*W; clear W

% rellena datos faltantes con interpolacion lineal
idn=find(isnan(sum(sh)));

%tic
%if any(isnan(sh(:)))
if ~isempty(idn)
    x=1:n;
    for i=idn; %i=1:p;
        t=sh(:,i);
        ai=find(~isnan(t));
        if length(ai)>1
            sh(:,i)=interp1(x(ai),t(ai), x)';
            %toc
        end
    end
    sh(isnan(sh))=0;
end

% Slepian tapers
  [psi] = dpss(n,nw,k);

% padding with zeros
  if npad==0; npad=2^(nextpow2(n)); end

% Despliega parametros de la estimacion espectral
  disp([' Slepian Tapers :        ' num2str(k) ])
  disp([' Bandwidth parameter :   ' num2str(nw) ])
  disp([' Padding :               ' num2str(npad) ])
  disp(['Performing ....' ])

  nf= npad/2;  % dimensiones de la mitad del espectro
  ddf=1./(npad*dt); % frecuencia de muestreo en npad intervalos
  fr  = ([1:nf]-1')*ddf; % vector de frecuencias

% Realiza la estimacion espectral
 nev=[];
 for i1=1:k
    tic
    psimat=repmat(psi(:,i1),1,p);   % Redimensionando tapers
    psimat=fft(psimat.*sh,npad);    % Estimativos spectrales con tapers
    psimat(nf+1:end, :)=[];         % Conservando mitad del espectro
    skn=['sk' num2str(i1)];         % Nombre matriz de almacenamiento
    eval([skn '= psimat;']);        % Matriz
    nev=[nev skn '(j1, :)''  '];    % matriz de evalaucion
    disp([' >Spectral stimated ' num2str(i1) ': ' num2str(toc) ' s.'])
 end
 clear psimat

% SVD de la matriz de espectros para diferentes frecuencias
  LFVS=fr*nan;
  tic
  for j1=1:nf
    eval(['M=[' nev '];'])
    [S] = svd(M, 0);
    % determine local fractional variance in the frequency fr(j1)
    %LFVS(j1) = S(1).^2/sum(S(2:end).^2); 
    LFVS(j1) = S(1).^2/sum(S(1:end).^2); 
  end
  vfbw=2*nw*(1/dt)*(1/n); % ancho de banda de la fecuencia espectral

% almacena los resultados en un archivo estructurado  
  LFV.name='Local Fractional Variance Spectrum'; 
  LFV.timeinterval=dt;
  LFV.bandwidth=nw;
  LFV.tapers=k;
  LFV.pading=npad;
  LFV.varfreqbandwidth=vfbw;
  LFV.spectrum=LFVS;
  LFV.specdomain=fr;

 disp([' >Local Fractional Variances : ' num2str(toc) ' s.' ])
 disp(['Done.' ])

