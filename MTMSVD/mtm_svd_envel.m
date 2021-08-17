function [env]=mtm_svd_envel(ff0,iif,fr,dt,ddf,n,k,psi,V,jsmoo,imode)
% Estima la funcion envolvente que se ajusta a un grupo de series de datos,
% a partir de la componente espectral (V) proveniende de la descomposicion
% MTM-SVD para una frecuencia fo. 
%
% Programa traducido, adecuado y optimizado en Matlab, a partir de codigos
% originales Mann&Park (1999)
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
% 2012.12.10  (ver. 1.0)

% save A
% clear all; load A;
% jsmoo=2;

ex=ones(1,n);

% La señal dominate de x y la envolvente son formalmente identicas para los
% modes seculares, donde df1=0;
df1=0; 
% df1=ff0/2

%--------------------------
%[n,c,s]=cossin(n,c,s,c0,s0,(df1),(dt));
c0=1; s0=0;
c(1)=c0;
s(1)=s0;
cs=cos(2.*pi.*df1.*dt);
sn=sin(2.*pi.*df1.*dt);
% cs=cos(2.*pi.*ff0.*dt);
% sn=sin(2.*pi.*ff0.*dt);
for  i=2:n;
    c(i)=c(i-1).*cs-s(i-1).*sn;
    s(i)=c(i-1).*sn+s(i-1).*cs;
end
%figure; plot(c);hold on; plot(s); title('kernels')
%--------------------------

% load in spectral envelopes for each mode.
%  multiply by 2.0./inc to obtain proper normalization;
%  2 from <cos(t).^2>=fix(1./2);
%  and inc from decimation of tapers in the data kernels;
%>>>  d=complex(ta(n1,:), tai(n1,:))*2; 

d=V(:,1).'; d=conj(d).*2; % doble del complejo conjugado de V

if iif==1
   d=V(:,1).'; d=conj(d); % complejo conjugado de V en la banda secular
end
%figure; plot(real(d) ); hold on;  plot(imag(d), 'r')

%  mult the tapers by sinusoid to get kernels;
%  note that we use the complex conjugate;
g=[];
for i0=1:k;
    g=[g ex(:).*complex(psi(:,i0).*c(:), -psi(:,i0).*s(:)) ];
end
%figure; plot(real(g) ); hold on; plot(imag(g), 'r')

% % --------------------------
% %[za,g]=premult(iboost,jsmoo);
 feps=1e-4;
 % dot prod with decaying envelope for fixed element;if we have boosted the
% kernels, the fixed element is constant 
c=1;
za=conj(sum(g*c));

%  then integrate the kernels jsmoo times
if(jsmoo>0)
    for ksmoo=1:jsmoo;
         g=flipud(cumsum(flipud(g))); % integracion de delante hacia atras
         
         %  little fudge factor to ensure penalty of first derivative;
         %g(1,:)=g(1,:)./feps; %OJO
    end
end
% figure; plot(g, 'r'); title(['kernel integrada ' num2str(jsmoo) ' veces'])
%  set(gca, 'ylim', [-30 30])

%--------------------------
% Descomposicion ortogonal de los tapers
%[g0,qrsave0, ipvt ] = zqrdc1(g,n,n,k,0,0); a=triu(g0,0); a(1:4, :);
%figure; plot(g0); title('descomposicion qr');set(gca, 'ylim', [-.08 .08])

[g1,qrsave1] = qr(g, 0); %qrsave=diag(qrsave)
%figure; plot(g1); title('descomposicion qr');set(gca, 'ylim', [-.08 .08])
%--------------------------

%--------------------------
%solve for the constant term;
%  first we calculate double backtransform of data vector (d)
%  and the triangular matrix r;      [X = R\(R'\B)]
dum = (conj(qrsave1)\(conj(qrsave1')\d.')).'; % dum=zbacktr1(k,n,g0,d,2)
%figure; plot(real(dum), 'b') ; hold on;  plot(imag(dum), 'r')

%  mult by za to get a scalar;
amp0=sum(conj(za).*dum); %! note: conjugated

%  next we calculate double backtransform of za;
%  and the triangular matrix r;      [X = R\(R'\B)]
dum = (conj(qrsave1)\(conj(qrsave1')\za.')).' ;%[dum]=zbacktr1(k,n,g0,za,2)

%  mult by za to get a scalar;
amp1=sum(conj(za).*dum); %! note: conjugated

% relacion de amplitudes
amp0=amp0./amp1;

sum1=sum(abs(d).^2);

d=d-za.*amp0;

sum2=sum(abs(d).^2);

rat=sum2./sum1;

%  solve forml m-tilde;
%  we backtransform once, then multiply by q; [X = U'\B]
env0=((conj((qrsave1)')\d.')).' ; %[env0]=zbacktr1(k,n,g0,d,1)  
%figure; plot(real(env0), 'b') ; hold on;  plot(imag(env0), 'r')
%--------------------------

%--------------------------
% recomponiendo envolvente
%zqrsl(g,n,n,k,qrsave,env,env,dum,dum,dum,dum,10000,info);
%env2=env0; env2(end+1:n)=0; [ env ] = zqrsl1(g0,n,n,k,qrsave0,env2,10000);
env = g1*(env0.');
%figure; plot(real(env), 'k') ; hold on;  plot(imag(env), 'r')
%--------------------------

%--------------------------
% Re-integrando envolvente
%
%[jsmoo,n,env,amp0,c]=postmult(jsmoo,n,env,amp0,c);
%function [jsmoo,npts,env,amp0,cs]=postmult(jsmoo,n,env,amp0,cs);
feps=1.0e-4;
if jsmoo>0
    %  post normalize to the envelope fluctuation;
    for ksmoo=1:jsmoo;
        %env(1)=env(1)./feps;
        env=cumsum(env);
    end
end
%  add the decay envelope back in;
env=env.'+amp0.*c;
npts2=n.*2;
% figure; plot(real(env), 'g') ; hold on;  plot(imag(env), 'r')
%--------------------------


