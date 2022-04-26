function GF=lowpassfilter(GFNS,dt,tau,k)
GF=zeros(1,k);
alpha=dt/(tau+dt);
GF(1)=alpha*GFNS(1);
for i=2:k
    GF(i)=alpha*GFNS(i)+(1-alpha)*GF(i-1);
end

end
