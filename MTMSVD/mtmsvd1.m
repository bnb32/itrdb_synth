%Written 02/06 by Toby Ault
%
%Last Modified: 20/02/06
%
%Matlab routine to perform MTM-SVD as described by Mann and Park, 1999 (1).
%This script expects a matrix whose rows are observations in time and whose
%columns are individual sites (i.e. gridpoints, sites of proxy data,
%etc...). It calculates and plots the Local Fractional Variance (lfv) 
%spectrum and stores the information in the variable lfv.
%
%Default values for K and p are 3 and 2, respectively (see Mann and Park, 
%1999). 
%
%References:
%
%(1) Mann, M.E., and Park, J., (1999) Oscillatory spatiotemporal signal
%detection, "Advances in Geophysics",41,1-131

clear
close all

% load yr_opt.mat
% tvalid = 200:401;
% X = V2.tsm(tvalid,:);

load synthdatn2.mat
X = data.tsm;

%Some generally useful variables for analysis later:
[N,M] =size(X); 
p = 2; %Time-Bandwith parameter for MTM
K =2*p-1; %Number of Data Tapers
nf = 4*N; %Number of frequencies to compute with Fourier Transform
uf = ceil(nf/2); %length of useful frequencies from Fourier transform
f = [0:nf]/(nf); %Frequencies of Fourier transform;
fr = 1./N; %Rayleigh Frequency
fbw =2*p*fr; %Frequency bandwidth resolution
f0bw =fbw/2; %Frequency bandwith from f0
fbwl = length(find(f <= fbw)); %Length of Bandwidth in f
f0bwl = length(find(f <= f0bw)); %Length of bandwith out from f0
per = 1./f; %Period

%First, normalize data and convert NaN's to Zeros:
for i = 1:M;
    q = find(~isnan(X(:,i)));
    xmean(i) = mean(X(q,i));
    X(q,i) = X(q,i)-xmean(i);
    ystd(i) = std(X(q,i));
    X(q,i) = X(q,i)./ystd(i);
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
    clear junk wmx Xtpr
    [junk, wmx] = meshgrid(1:M, w(:,i));
    Xtpr = wmx.*X;
    dft(:,:,i) = fft(squeeze(Xtpr),nf);
end

%Perfrom SVD for every frequency, f, on the matrix Y(k,f) whose rows are 
%spatial locations and whos columns are the K Fourier transformed 
%tapered data.

for i = 1:uf;
    clear Yk Utemp Vtemp Ltemp
    for j = 1:K
        Ykf(:,j) = squeeze(dft(i,:,j))';
    end    
    [Utemp, Ltemp, Vtemp] = svd(Ykf);
    U(i,:) = Utemp(:,1) ;
    V(i,:) = Vtemp(:,1);
    sv = diag(Ltemp).^2;
    lfv(i) = sv(1)/sum(sv);
    L(i,:) = sv;
end
 
plot(lfv(1:uf));
hold on
% plot(1/50, [0:.01:1])
%  plot(1/30, [0:.01:1])
%  plot(1/10, [0:.01:1])

%sigmtmsvd1