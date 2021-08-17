function Ykf=parsqueeze(dft,M,K,i);
    Yfk=nan(M,K);
    for j = 1:K
        a=squeeze(dft(i,:,j))';
        Ykf(:,j) = parmat(M,[1:M],a(:));%squeeze(dft(i,:,j))';
    end
end    
