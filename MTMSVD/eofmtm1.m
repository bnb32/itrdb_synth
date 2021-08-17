
1/f(f0)

for i = 1:k
    A(i,:) = (L(f0,i)*conj(V(f0,i)).*w(:,i)')./c(i);
end
A1 = sum(A);
Arms = std(A1);

if f0 >= fbwlength 
    yf0 = 2;
elseif f0 < fbwlength
    yf0 = 1;
end

eofm = yf0*ystd.*(U(f0,:))*Arms; %Spatial pattern
sgnlm = yf0*real(ystd(m).*U(f0,m).*(A1).*exp(-1i*2*pi*f(f0)*[1:N]));