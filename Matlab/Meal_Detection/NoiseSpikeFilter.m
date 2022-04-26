function GFNS=NoiseSpikeFilter(Gm,dG)
%%%%%%%%%%%%%%%
%Function to calculate the Noise-Spike filter
%% input
%[Vector] Gm Measure ment (Length: k)
%[int] dG Maximum allowable ROC
%% Output
%[Vector] GFNS with length k Noise-spikefilter has been applied
k=length(Gm);
GFNS=zeros(1,k);

GFNS(1)=Gm(1);
for i=2:k
    if abs(Gm(i)-GFNS(i-1))<=dG
        GFNS(i)=Gm(i);
        
    elseif (GFNS(i-1)-Gm(i))>dG
        GFNS(i)=GFNS(i-1)-dG;
        
    elseif (Gm(i)-GFNS(i-1))>dG
        
        GFNS(i)=GFNS(i-1)+dG;
    end
end
end