function dGF=LagInterpol(GF,t)
%%NOT IN USE
 k=length(GF);
 dGF=zeros(1,k);
 for i=1:k
    dGF(i)=(t(i)-t(i-1))/((t(i-2)-t(i-1))*(t(i-2)-t(k)))*GF(i-2)...
        +(t(i)-t(i-2))/((t(i-1)-t(i-2))*(t(i-1)-t(i)))*GF(i-1)...
        +(2*t(i)-t(i-2)-t(i-1))/((t(i)-t(i-1))*(t(i)-t(i-2)));
end
end