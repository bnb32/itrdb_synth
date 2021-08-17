function [lfv,stats] = parsvd(X,varargin)

addpath('../');
%
%Written 02/06 by Toby Ault
%Updated by Brandon Benton
%
%Last Modified: 06/07/2019
%
%Usage:
%stat = mtmsvdf1(X,varargin)
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

X_tmp=nan(N,M);
%ystd=nan(1,M);
%xmean=nan(1,M);

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

dt=1/12;
p = 2; %Time-Bandwith parameter for MTM
K =2*p-1; %Number of Data Tapers
uf = ceil(nf/2); %length of useful frequencies from Fourier transform
f = [0:nf]/(dt*nf); %Frequencies of Fourier transform;
fr = 1./(dt*N); %Rayleigh Frequency
fbw =2*p*fr; %Frequency bandwidth resolution
f0bw =fbw/2; %Frequency bandwith from f0
fbwl = length(find(f <= fbw)); %Length of Bandwidth in f
f0bwl = length(find(f <= f0bw)); %Length of bandwith out from f0
frbwl = length(find(f <= fr));

%Normalize data and convert NaN's to Zeros:

if norm == 'Y'
    %disp(['Normalizing data for ' num2str(M) ' sites: ']);
    parfor i = 1:M;
        
	q = find(~isnan(X(:,i)));
        xmean(i) = mean(X(q,i));
	a=X(q,i)-xmean(i);
	ystd(i)=std(a);
	b=a./ystd(i);
	X_tmp(:,i)=parmat(N,q,b);
    end

elseif norm == 'N'
    xmean=0;
    ystd=1;
    Xtmp=X(:,:);
else
    error('Normalization command not recognized')
end

%Convert NaN's to Zero. This works because the Discrete Fourier Transform
%assumes that the N observations are sampled from a continous process. 
%Zeros at the ends or between observations will have no impact on the
%frequencies derived.
X_tmp(find(isnan(X_tmp))) = 0;

%Compute Discrete Prolate Slepian Sequences:
[w, c] = dpss(N,p,K);

%for k=2:K
%    w(:,k)=(-1)*w(:,k);
%end    

%Build nf x M x K matrix of Fourier transforms for each tapered dataset 

%disp(['Building fourier transform matrix of ' num2str(K) ' tapered datasets: ']);

%if isempty(gcp('nocreate'));
%    ncores=feature('numcores');
%    parpool(ncores);
%end

parfor k = 1:K
    [~, wmx] = meshgrid(1:M, w(:,k));
    Xtpr = wmx.*X_tmp;
    dft(:,:,k) = fft(Xtpr,nf);
    %for i=1:nf
    %    emat=exp(-1i*2*pi*f(i)*[1:N]);
    %    for m=1:M
    %        dft(i,m,k)=dot(Xtpr(:,m),emat(:));
    %	end
    %end	
end

%Perfrom SVD for every frequency, f, on the matrix Y(k,f) whose rows are 
%spatial locations and whos columns are the K Fourier transformed 
%tapered data.

%disp(['Performing SVD for ' num2str(uf) ' frequencies: ']);

for i = 1:uf;
    for j = 1:K;
        Ykf(i,:,j)=squeeze(dft(i,:,j))';
    end
end    

parfor i = 1:uf;
    Ykf_tmp=squeeze(Ykf(i,:,:));
    [Utemp, Ltemp, Vtemp] = svdecon(Ykf_tmp); 
    U(i,:,:)=Utemp(:,:);
    V(i,:,:)=Vtemp(:,:);
    %U(i,:) = Utemp(:,1);
    %V(i,:) = Vtemp(:,1);
    sv = diag(Ltemp);
    ev = diag(Ltemp).^2;
    svals(i)=sv(1);
    lfv(i) = ev(1)/(sum(ev));
    L(i,:) = sv;

%% is this doing what i think?
%% why use Utemp(:,1) instead of Utemp(:,:), etc?

    %U(i,:,:)=Utemp(:,:);
    %V(i,:,:)=Vtemp(:,:);
    %L(i,:,:)=Ltemp(:,:);


end


%disp(['Finished SVD for ' num2str(uf) ' frequencies: ']);
 
%sigmtmsvd1
stats.lfv = lfv;
stats.sv = svals;
stats.L = L;
stats.U = U;
stats.V = V;
stats.p = p; %Time-Bandwith parameter for MTM
stats.K =K; %Number of Data Tapers
stats.nf = nf; %Number of frequencies to compute with Fourier Transform
stats.uf = uf;%f(1:uf); %length of useful frequencies from Fourier transform
stats.f = f; %Frequencies of Fourier transform;
stats.fr = fr; %Rayleigh Frequency
stats.fbw = fbw; %Frequency bandwidth resolution
stats.f0bw = f0bw; %Frequency bandwith from f0
stats.fbwl = fbwl; %Length of Bandwidth in f
stats.f0bwl =f0bwl; %Length of bandwith out from f0
stats.frbwl = frbwl;
stats.meanm = xmean; %Mean of each site
stats.stdm = ystd;   %Standard deviation of each site
stats.X = X_tmp;
stats.Y = Ykf;
