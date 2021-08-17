function [stats]=get_signal(stats);

[N,M]=size(stats.X);

uf=stats.uf;
K=stats.K;
U=stats.U;
V=stats.V;
L=stats.L;
p=stats.p;
ystd=stats.stdm;
f=stats.f;
Y=stats.Y;


U=U(:,:,1);%squeeze(mean(U,3));
V=V(:,:,1);%squeeze(mean(V,3));

[w,c]=dpss(N,p,K);

%for k=2:K
%    w(:,k)=(-1)*w(:,k);
%end    

A=[];
Asum=[];
eof=[];
signl=[];
ts=[];

%dt=1/12.0;

%for k=1:K
%    cinv(k)=1.0/c(k);
%end    

wmat=nan(N,K);
lvmat=nan(uf,K);
B=nan(N,M,K);
Asum=nan(uf,N);


if isempty(gcp('nocreate'));
    ncores=feature('numcores');
    parpool(ncores);
end


%parfor n=1:N
%    wmat(n,:)=w(n,:).*cinv(:)';
%end

%parfor i=1:uf
%    lvmat(i,:)=L(i,:).*conj(V(i,:));
%end    

%parfor k=1:K
    %for i=1:uf
    %    for n=1:N
%            A(:,:,k)=lvmat(:,k)*wmat(:,k)';
%	end
%    end	
%end

%Asum=sum(A,3);

%parfor m=1:M;
%    LinvY(:,m,:)=Linv(:,:).*squeeze(Y(:,m,:));
%end    
countDir='/home/bnb32/projects/itrdb_synth/MTMSVD/countDir';
psimple_progress(countDir,0,uf);


parfor i=1:uf
    psimple_progress(countDir,i,uf)
    [B,nmat]=Bfill(L,Y,w,N,M,K,f(i),i);
    %for n=1:N;
    %   for m=1:M;
    %       for k=1:K;
    %           B(n,m,k)=1.0/(L(i,k))*Y(i,m,k)*w(n,k);
    %       end
    %   end
    %end


    temp=sum(B,2);
    temp=squeeze(sum(temp,3));
    Asum(i,:)=temp;
    Arms(i)=std(Asum(i,:));
    %max(Asum(i,:));%std(Asum(i,:));
    mmat=ystd.*U(i,:);
    eof(i,:)=mmat*Arms(i);
    exps(i,:)=exp(-1i*2*pi*f(i)*[1:N]);
    %nmat=Asum(i,:).*exps(i,:);
    %signl(i,:,:)=real(nmat'*mmat);
    ts(i,:,:)=real(nmat);
end

stats.ts=ts;
stats.eof=eof;
