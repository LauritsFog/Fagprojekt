function GRID=GridAlgo(Gm,dG,dt,tauF,t)
%%%%%%%%%%%%%%%
%Function to calculate the Grid algorithm
%% input
%[Vector] Gm Measure ment (Length: k)
%[int] k length of vectors
%[int] dG Maximum allowable ROC
%% Output
%[Vector] GRID with length k either 0 or 1 at the i'th entry.
%% Inisializing
k=length(Gm);

dGF=zeros(1,k);
%% preprocessing section
GFNS=NoiseSpikeFilter(Gm,dG);

%% Low pass filter

GF=lowpassfilter(GFNS,dt,tauF,k);

%% Lagrangian interpolation
%PROBABLY WRONG NEEDS SOME INITIAL VALUES
for i=1:k
    T1=(t(i)-t(i-1))/((t(i-2)-t(i-1))*(t(i-2)-t(i)));
    T2=(t(i)-t(i-2))/((t(i-1)-t(i-2))*(t(i-1)-t(i)));
    T3=(2*t(i)-t(i-2)-t(i-1))/((t(i)-t(i-1))*(t(i)-t(i-2)));
    dGF(i)=T1*GF(i-2)...
        +T2*GF(i-1)...
        +T3*GF(i);
end

%% GRID
GRID=zeros(1,k);
for i=1:k
    if GF(i)>Gmin && ((dGF(end-2:end)>Gmin3) || (dGF(end-1:end)>Gmin2))
        GRID(i)=1;
    else
        GRID(i)=0;
    end
end



end