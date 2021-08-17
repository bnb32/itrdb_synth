function [R,RC,RP,U,S,V]=mtm_svd_bandrecon(sh,nw,k,dt,freq,imode,npad,lan,vw)
%function [R,RC,RP,U,S,V]=mtm_svd_bandrecon(sh,nw,k,dt,freq,imode,npad,lan,vw)%
%
%    Perform spatial & temporal mode reconstructions at a particular
%    frequency for a peak detected by MTM-SVD method
%
%
% Inputs:
%
%      sh       Time series matrix [rows=time columns=space]
%      nw       Time-frequency bandwidth parameter [2]
%      k        Number of Slepian Tapers. Just the k=2*nw-1 are resistant
%               to spectral leakeage [3]. 
%      dt       Time interval (cycles year-1). Instead of the time interval
%               parameter, a regular date vector in Matlab format (as the
%               produced by datenum.m) can be entered.
%   freq        Frequencies to reconstruc (cycles year-1) [1xF]
%   imode       Number of the mode to reconstruct [1]. 
%   npad        Set frequency range for DFT, frequency spacing, zero
%               padding, etc. npad must be a power of 2 greater than or
%               equal to the time series length (nscan). Larger values
%               provide a more interpolated frequency axis. If npad is left
%               equal to 0, npad will set like 2^(nextpow2(nscan)+1).
%   lan         Inversion option for time domain:
%               lan==0, "Minimum Norm" : Minimized the overall size of the
%               envelope of the oscillation. Generally is the most faithful
%               constraint, but may artificially minimize amplitude of the
%               ends of the time interval.  
%               lan==1, "Minimum Slope" : Minimized the average slope of
%               the envelope of the oscillation. Particularly useful for
%               reconstructing secular trends in the series which do not
%               have the periodic properties of the non-secular
%               oscillations.    
%               lan==2, "Minimum Roughness" : minimize the average 2nd
%               derivative of the envelope of the oscillation. Useful if
%               the amplitude of the envelope is believed to be changing
%               abruptly near the ends of the time interval.
%               lan>2, "Minimum Error" : minimize the mean squared error
%               between raw data and the reconstruction over all possible
%               linear combinations of the alternative apriori constraints
%               described above. Methods available:
%               lan==3, min Norm/min Roughness 
%               lan==4, min Norm/min Roughness/min Slope
%   vw          Weights on dataseries: the user may decide to set e.g., 
%               areal(vw=cos(latitude) [space]) weights on the constituent
%               dataseries. Otherwise uniform weights are set if vw=0
%
% Outputs:  
%
%   U, S, V     Martrices of complex spatial orthogonal functions [M x k],
%               singular values [k x 1] and principal components (or
%               modulations) [k x k]  respectively.
%   R           Reconstructed matrix [MxN]
%   RC          Reconstructed cycle [Nx360xF]
%
%   RP          Structure array containig:
%   RP.recfrequencies       (Reconstructed frequencies)
%   
%   RP.timeinterval         Time interval (cycles/year)
%   RP.bandwidth            Bandwidth parameter used
%   RP.tapers               Slepian Tapers used
%   RP.pading               Pading in the spectrum domain
%   RP.amplitude            (Peak amplitude variation at each site [MxF])
%   RP.phase                (relative to first position [MxF])
%   RP.fractionalvar        (Total fractional variance for each mode [kx1])
%   RP.partialvar           (Partial fractional variance for each mode [kx1])
%   RP.varexplained         (Percent of raw variance explained for the
%                            reconstructed mode on each position [Mx1])
%   RP.totvarexplained      (Total variance explained)
%   RP.lambda               (Linear combination of the envelopes functions
%                            used in the reconstruction [1xk])
%   RP.temporalmean         (mean field)
%   RP.temporalstd          (standard deviation field)
%   RP.eof                  (spatial eofs)
%
%
%--------------------------------------------------------------------------
% A description of the theoretical basis of the MTM-SVD toolbox and some
% implementation details can be found in:
%
%  Correa-Ramirez, M. & S. Hormazabal, 2012. "MultiTaper Method-Singular
%  Value Decomposition (MTM-SVD): Variabilidad espacio–frecuencia de las
%  fluctuaciones del nivel del mar en el Pacífico Suroriental". Lat. Am. J.
%  Aquat. Res., 40[4]: 1039-1060 
%
% Citation of this article would be appreciated if you use the toolbox.
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
% 2014.09.10  (ver. 2.0)

  [n, p]=size(sh); 

% measurement interval
  if length(dt)>1
    lt=(dt(end)-dt(1))/365.25;    % series length in years
    dt=lt/length(dt);  
  end

% standardize data and select envelope method
  vm = nanmean(sh);
  vs = nanstd(sh);
  sh = sh-repmat(vm,n,1);
  sh = sh./repmat(vs,n,1);

  if lan==0; lambda=[1 0 0]; end
  if lan==1; lambda=[0 0 1]; end
  if lan==2; lambda=[0 1 0]; end
  if lan==3; lambda=[1 1 0]; end
  if lan==4; lambda=[1 1 1]; end

% spatial location weights
if vw==0; vw=ones(1,p); end

% multiply data by spatial weights
  W = repmat(vw(:)',n,1);
  %whos W sh
  sh = sh.*W; clear W

% Define padding for DFT
  if npad == 0; npad = 2^(nextpow2(n)); end

  ddf=1./(npad*dt);      % frequency interval 
  nf= npad/2;            % length of usable frequencies 
  fr  = ([1:nf]-1')*ddf; % frequencies vector

% reconstruction parameters
  disp(['MTM-SVD BAND RECONSTRUCTION   '])
  disp([' Slepian Tapers :        ' num2str(k) ])
  disp([' Bandwidth parameter :   ' num2str(nw) ])
  disp([' Padding :               ' num2str(npad) ])
  disp(['Performing ....' ])

% fill in missing data with linear interpolation 
idn=find(isnan(sum(sh)));   % missing data indices
if ~isempty(idn)
    x=1:n;
    for i=1:length(idn);
        t=sh(:,idn(i));
        ai=find(~isnan(t));
        if length(ai)>1
            sh(:,idn(i))=interp1(x(ai),t(ai), x)';
        end
    end
end
sh(isnan(sh))=0;

% construct slepian tapers 
[psi] = dpss(n,nw, k);

% get the spectral estimates for the slepian tapers
  nev=[];
  for i1=1:k
    tic
    psimat=repmat(psi(:,i1),1,p);   % taper matrix
    psimat=fft(psimat.*sh,npad);    % spectral estimation for tapers
    psimat(nf+1:end, :)=[];         % keep half of the spectrum
    skn=['sk' num2str(i1)];         % storage array name
    eval([skn '= psimat;']);        % save the spectrum matrix
    nev=[nev skn '(iif, :)''  '];   % name of evlaluation matrix
    %disp([' > Spectrum ' num2str(i1) ' estimated in ' num2str(toc) ' seconds' ])
  end
  %whos psimat fr
  clear psimat skn

% Define the output matrices
  S = ones(k, length(freq))*nan; fracvar = S; partvar = S;
  R1 = sh'*0; R2 = R1; R3 = R1;
  
  %disp(['Reconstructed frequencies : ' ]);
  D = repmat(vs,n,1)';     % standard deviation matrix
  W = repmat(vw(:)',n,1)'; % spatial weight matrix

  for i1=1:length(freq);
      % index of the closest frequencies
      iif = near(fr, freq(i1)); if iif==0; iif=1; end
      
      % nearest frequency
      ff0 = fr(iif);
      %disp(['(' num2str(i1) ') ' num2str(ff0) ' cycles/yr | ' num2str(1/ff0) 'yr']);
      
      % build matrix used in SVD 
      eval(['M=[' nev '];']);
      [U(:,:,i1),S0,V(:,:,i1)] = svd(M, 0);
      S(:,i1) = diag(S0);
     
      % calculate lfv for each mode
      fracvar(:,i1)=S(:,i1).^2./sum(S(:,i1).^2);
      
      % calculate partial lfv for each mode
      for j=1:k
          partvar(j,i1)=(S(j).^2)./sum(S(j:end).^2);
      end
      
      % calculate envelope based on selected constraints

      % minimize norm
      [env1]=mtm_svd_envel(ff0,iif,fr,dt,ddf,n,k,psi,V(:,:,i1),0,imode);
      
      % minimize 2nd derivative
      [env2]=mtm_svd_envel(ff0,iif,fr,dt,ddf,n,k,psi,V(:,:,i1),2,imode);

      % minimize 1st derivative
      [env3]=mtm_svd_envel(ff0,iif,fr,dt,ddf,n,k,psi,V(:,:,i1),1,imode);

%      figure; plot(real(env1)); hold on; plot(real(env2)); plot(real(env3))
%      V
%      hold on; 
%      plot(real(cumsum(psi*V(:,1))), 'k')
     
     % calculate sinusoids for reconstruction over time
      cs=1;
      sn=0;
      c=cos(2*pi*ff0*dt);
      s=sin(2*pi*ff0*dt);
      
      for i2=2:n;
          cs(i2)=cs(i2-1)*c-sn(i2-1)*s;
          sn(i2)=cs(i2-1)*s+sn(i2-1)*c;
      end
      CS = complex(cs.', sn.'); CS=conj(CS);
      %CS = complex(cs, sn); CS=conj(CS)*-1;
 
%       figure; plot(cs)
%       

      % Maximum amplitude of the envelopes (complex) attained over record
      A=[(env1) (env2) (env3)];[a,b]=max(abs(A)); envmax(i1)=real(A(b));
       %env3= envmax;   env2= envmax;      env1= envmax;

      %if lan==0;
      %    envmax(i1)=std(env1);
      %elseif lan==1;
      %    envmax(i1)=std(env3);
      %else;
      %    envmax(i1)=std(env2);
      %end	  
       
      %-----------------------------------------------------------------
      % RECONTRUCTIONS
      
      % R1 "Minimum Norm"
      % Minimized the overall size of the envelope of the oscillation.
      % Generally is the most faithful constraint, but may artificially
      % minimize amplitude of the ends of the time interval.
      
      % R3 "Minimum Slope" :
      % Minimized the average slope of the envelope of the oscillation.
      % Particularly useful for reconstructing secular trends in the series
      % which do not have the periodic properties of the non-secular
      % oscillations.
      
      % R2 "Minimum Roughness" :
      % Minimized the average of the 2nd derivative of the envelope of the
      % oscillation. Useful if the amplitude of the envelope is believed to
      % been changing abruptly near the ends of the time interval.
      
      % Do time reconstructions for indicated mode for all sites indicated
      % Take complex conjugate of carrier frequency sinusoid and do the
      %  take Re{ std  A(t) U_l(x) exp(iw_p t) } to reconstruct oscillation
      %  at grid point
      %
      % A=the matrix produc U x the root mean squared amplitud of the
      % demodulating enveloped A1(fo,t).
      
      
      R1 = real( D.*S(imode, i1).*(U(:,imode, i1)*(CS.*env1.').')./W );
      
      %R1 = R1 + real( S(imode, i1).*(U(:,imode, i1)*(CS.*env1.').') );
     
      R2 = R2 + real( D.*S(imode, i1).*(U(:,imode, i1)*(CS.*env2.').')./W );
      R3 = R3 + real( D.*S(imode, i1).*(U(:,imode, i1)*(CS.*env3.').')./W );
      
  end
  
  clear W D M sk*
  
  % Case min norm/min roughness/min slope
  if lan == 4
      disp(['Error minimization of the Reconstruction : ' ]);
    
    % create vector to save time
    shr=sh'.*repmat(vs,n,1)';
    
    % linear input vector
    L=0:.1:1;
    L=nchoosek(L,2);
    L(:,3)=1-sum(L,2); L(L(:,3)<0, :)=[];
    
    % possible permutations
    plan= perms(1:3);
    
    % linear combinations
    lcom=[];
    for i1=1:size(plan,1);
        lcom=[lcom; L(:,plan(i1,:))];
    end       
    
    % remove repeated combinations
    lcom = unique(lcom, 'rows');
    
    % output matrix
    SER = lcom(:,1)*nan;
    for i1=1:length(lcom)
        tic
        % reconstruction
        R=[];
        R=lcom(i1,1).*R1 + lcom(i1,2).*R2 + lcom(i1,3).*R3;        
        R=(shr-R).^2; % sum of squares of the reconstruction error
	SER(i1)=nansum(R(:));
        if i1==1
            tt=(toc*(length(lcom)-i1))./60;
            %disp(['     Time to complete : ' sprintf('%4.1f', tt) ' min' ])
            %disp(['     processing ... ']); pause(.1)
        end
    end
    clear shr
    [a,b] = min(SER); % select combination with lowest error
    lambda = lcom(b,:);
    disp(['  Done.' ])
end

% Case min norm/min roughness 
if lan ==3 
    disp(['Error minimization of the Reconstruction :' ]);
    shr=sh'.*repmat(vs,n,1)';
    
    % possible permutations
    lcom=[0:.01:1]';
    lcom(:,2)=1-lcom;
    
    SER=lcom(:,1)*nan; 
    for i=1:length(lcom)
        tic
        % reconstruction
        R=lcom(i,1).*R1 + lcom(i,2).*R2;
        R=(R-shr).^2; % sum of squares of reconstruction error
        SER(i)=nansum(R(:));
        if i==1
            tt=(toc*(length(lcom)-i))./60;
            %disp(['     Time to complete : ' sprintf('%4.1f', tt) ' min' ])
            %disp(['     processing ... ']); pause(.1)
        end
    end
    clear shr
    [a,b]=min(SER); % select combination with lowest error
    lambda=[lcom(b,:) 0];
    %disp(['  Done.' ])
end

% perform the final reconstruction with the smallest error coefficients
  R = lambda(1).*R1 + lambda(2).*R2 + lambda(3).*R3;  clear R1 R2 R3
  disp(['Inversion option used in the Reconstruction :' ])
  disp([' lambda = [ ' num2str(lambda) ' ]' ])


% CYCLE RECONSTRUCTION
%
% Determine spatial patterns, that corresponding to peak amplitude of the
% oscillation over the duration of record (note that this is 1/2
% peak-to-peak)     

% Peak amplitude variation at each site
% amp = Singular value * max|time envelope| * value of spatial eof 
%               * standard deviation_i / weight_i
 
% phase(amplitude(i)) = temporal phase relative to maximum amplitude at a
% reference site 

% Associated temporal lag = (phase/360)*1/freq
%disp(['Canonical cycle :' ]);
 snph = [0:1:360]'; % snapshop phases to be evaluated
 RC = nan(p,length(snph),length(freq)); %prelocating output matix
 %W = repmat(vw(:)',n,1)'; % weigths

for i1 = 1:length(freq) % canonical cycle for each frequency
    % Reference phase (First position)
    angle1 = 180.0*atan2(imag(U(1,imode, i1)), real(U(1,imode,i1)))/pi;
    
    % Phase (relative to first position)
    a=(180.*atan2(imag(U(:,imode,i1)),real(U(:,imode,i1)))/pi);%-angle1;
    %a(a<0)=a(a<0)+360;
    Rph(:,i1)=a;
    
    % Associated temporal lag (years)
    Rlag(:,i1)=Rph(:,i1)./(360.*freq(i1));
    
    % Amplitude
    amp=abs(U(:,imode,i1));
    
    % Peak amplitude variation at each site
    Ramp(:,i1)=S(imode, i1).*amp.*envmax(i1).*vs.'./vw(:);

    %spatial eofs
    %eof(:,i1)=Ramp(:,i1).*exp(1i*pi/180.0*Rph(:,i1));
    %S(imode,i1).*U(:,imode,i1).*envmax(i1).*vs.'./vw(:);     
    %-------------------------------------------------------
    % Cycle reconstructed
    for i=1:p
        tic
        Cph=snph-(Rph(i,i1));%+angle1);
        RC(i,:,i1)=Ramp(i,i1).*cos(Cph.*pi./180);        

        if i==1; tt=(toc*(p-i))./60; 
	    %disp(['     Time to complete : ' sprintf('%4.1f', tt) ' min' ]); 
	end
    end
end

%spatial eofs
[nx,ny]=size(Ramp);
for i=1:nx;
    for j=1:ny;
        eof(i,j)=Ramp(i,j)*exp(1i*pi/180.0*Rph(i,j));
    end
end
eof=mean(eof,2);

%eof=mean(Ramp,2).*exp(1i*pi/180.0*mean(Rph,2));


%-------------------------------------------------------

% reincorporate standard deviations and remove weights
sh=sh'.*repmat(vs,n,1)'./repmat(vw(:)',n,1)'; 

% % percent of raw variance explained
%  vexp=nan(1,p);
% for i=1:p
%     tic
%     b=corrcoef(R(i,:)', sh(i, :)).^2*100;
%     vexp(i)=b(2);            
% 
%     if i==1; tt=(toc*(p-i))./60; disp(['     Time to complete : ' sprintf('%4.1f', tt) ' min' ]); end
% end
% %toc
% a=find(~isnan(sh+R));
% a=corrcoef(R(a), sh(a)).^2*100;
% totvarexp=a(2);

%percent of raw variance explained
vsr=var(R');
vexp=vsr./(vs.^2)*100;
totvarexp=nansum(vsr)./nansum(vs.^2)*100;
%disp(['Variance explained in the Reconstruction : ' sprintf('%4.2f', totvarexp) '%']);

fray = 1/(n*dt);                    % rayleigh frequency
fbw=2*nw*fray;                      % rayleigh bandwidth
%whos

% output matrices
RP.name='MTM-SVD spatial and temporal frequency reconstruction'; 
RP.frequenciesrec=freq;
RP.timeinterval=dt;
RP.bandwidth=nw;
RP.tapers=k;
RP.pading=npad;
RP.varfreqbandwidth=fbw;
RP.amplitude=Ramp; clear Ramp
RP.phase=Rph;clear Rph
RP.lag=Rlag; clear Rlag 
RP.fractionalvar=fracvar;
RP.partialvar=partvar;
RP.varexplained=vexp; clear varexp
RP.totvarexplained=totvarexp;
RP.lambda=lambda;
RP.temporalmean=vm; clear vm
RP.temporalstd=vs; clear vs
RP.eof=eof; clear eof

disp(['Done.' ]) 
