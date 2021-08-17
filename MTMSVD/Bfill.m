function [B,nmat]=Bfill(L,Y,w,N,M,K,f,i);
    for n=1:N;
        for m=1:M;
            for k=1:K;
	        B(n,m,k)=1.0/(L(i,k))*Y(i,m,k)*w(n,k);
	    end
        end
    end

    temp=sum(B,3);
    for n=1:N;
        for m=1:M;
	    nmat(n,m)=temp(n,m,1)*exp(-1i*2*pi*f*n);
	end
    end	
