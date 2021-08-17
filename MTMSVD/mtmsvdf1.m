function [lfv,stat] = mtmsvdf1(X,varargin)
%
%Written 02/06 by Toby Ault
%
%Last Modified: 20/02/06
%
%Usage:
%[lfv,stat] = mtmsvdf1(X,varargin)
%
%Matlab function to perform MTM-SVD as described by Mann and Park, 1999 (1).
%This script expects a matrix whose rows are observations in time and whose
%columns are individual sites (i.e. gridpoints, sites of proxy data,
%etc...). 
%
%It returns the Local Fractional Variance (lfv) spectrum, along with
%spatial EOFs, "Principal Modulations," and singular values for each
%frequency.
%
%Default values for K and p are 3 and 2, respectively (see Mann and Park, 
%1999). 
%
%References:
%
%(1) Mann, M.E., and Park, J., (1999) Oscillatory spatiotemporal signal
%detection, "Advances in Geophysics",41,1-131

%Some generally useful variables for analysis later:
[N,M] =size(X); 
if nargin ==1
    nf =N*4;
    norm = 'Y';
elseif nargin ==2  
    nf = varargin{1};
    norm = 'Y';
elseif nargin ==3
    if isempty(varargin{1})
        nf = N*4;
    else
        nf = varargin{1};
    end
    norm = varargin{2};
else
    error('Invalid number of input arguments')
end

p = 2; %Time-Bandwith parameter for MTM
K =2*p-1; %Number of Data Tapers
uf = ceil(nf/2); %length of useful frequencies from Fourier transform
f = [0:nf]/(nf); %Frequencies of Fourier transform;
fr = 1./N; %Rayleigh Frequency
fbw =2*p*fr; %Frequency bandwidth resolution
f0bw =fbw/2; %Frequency bandwith from f0
fbwl = length(find(f <= fbw)); %Length of Bandwidth in f
f0bwl = length(find(f <= f0bw)); %Length of bandwith out from f0

%Normalize data and convert NaN's to Zeros:
if norm == 'Y'
    for i = 1:M;
        q = find(~isnan(X(:,i)));
        xmean(i) = mean(X(q,i));
        X(q,i) = X(q,i)-xmean(i);
        ystd(i) = std(X(q,i));
        X(q,i) = X(q,i)./ystd(i);
    end
elseif norm == 'N'
    X(find(isnan(X))) = 0;
    xmean =0;
    ystd =1;
else
    error('Normalization command not recognized')
end

%Convert NaN's to Zero. This works because the Discrete Fourier Transform
%assumes that the N observations are sampled from a continous process. 
%Zeros at the ends or between observations will have no impact on the
%frequencies derived.
X(find(isnan(X))) = 0;

%Compute Discrete Prolate Slepian Sequences:
[w, c] = dpss(N,p,K);

%Build nf x M x K matrix of Fourier transforms for each tapered dataset 
for i = 1:K
    clear wmx Xtpr
    [~, wmx] = meshgrid(1:M, w(:,i));
    Xtpr = wmx.*X;
    dft(:,:,i) = fft(squeeze(Xtpr),nf);
end

%Perfrom SVD for every frequency, f, on the matrix Y(k,f) whose rows are 
%spatial locations and whos columns are the K Fourier transformed 
%tapered data.

countDir='/home/bnb32/projects/itrdb_synth/MTMSVD/countDir';
psimple_progress(countDir,0,uf);

for i = 1:uf;
    clear Yk Utemp Vtemp Ltemp
    for j = 1:K
        Ykf(:,j) = squeeze(dft(i,:,j))';
        %[~, wmx] = meshgrid(1:M, w(:,j));
	%Akf(:,j) = wmx.*Ykf(:,j);
    end    
    %Ykf=squeeze(dft(i,:,:));

    psimple_progress(countDir,i,uf)
    [Utemp, Ltemp, Vtemp] = svdecon(Ykf);
    U(i,:) = Utemp(:,1) ;
    V(i,:) = Vtemp(:,1);
    sv = diag(Ltemp);
    ev = diag(Ltemp).^2;
    svals(i)=sv(1);
    lfv(i) = ev(1)/sum(ev);
    L(i,:) = sv;
end


A=[];
Asum=[];
eof=[];
signl=[];
ts=[];

%dt=1/12.0;

for i=1:uf
    for n=1:N
        for k=1:K
            %A(i,n,k)=sqrt(L(i,k))*conj(V(i,k))*w(n,k)./c(k);
            A(i,n,k)=(L(i,k))*conj(V(i,k))*w(n,k)./c(k);
        end    
        Asum(i,n)=sum(A(i,n,:));
	ts(i,n)=real(Asum(i,n)*exp(-1i*2*pi*f(i)*n));
    end
    Arms(i)=std(Asum(i,:));
    eof(i,:)=ystd.*(U(i,:))*Arms(i);
end

%sigmtmsvd1
stat.lfv = lfv;
stat.L = L;
stat.U = U;
stat.V = V;
stat.p = p; %Time-Bandwith parameter for MTM
stat.K =K; %Number of Data Tapers
stat.nf = nf; %Number of frequencies to compute with Fourier Transform
%stat.uf = f(1:uf); %length of useful frequencies from Fourier transform
stat.uf = uf; %length of useful frequencies from Fourier transform
stat.f = f; %Frequencies of Fourier transform;
stat.fr = fr; %Rayleigh Frequency
stat.fbw = fbw; %Frequency bandwidth resolution
stat.f0bw = f0bw; %Frequency bandwith from f0
stat.fbwl = fbwl; %Length of Bandwidth in f
stat.f0bwl =f0bwl; %Length of bandwith out from f0
stat.meanm = xmean; %Mean of each site
stat.stdm = ystd;   %Standard deviation of each site
stat.X=X;
stat.Xhat=eof;
stat.ts=ts;
