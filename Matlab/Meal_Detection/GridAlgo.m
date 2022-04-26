function [GF,dGF,GRID]=GridAlgo(Gm,dG,dt,~,t)
%%%%%%%%%%%%%%%
%Function to calculate the Grid algorithm
%% input
%[Vector] Gm Measure ment (Length: k)
%[int] k length of vectors
%[int] dG Maximum allowable ROC
%% Output
%[Vector] GRID with length k either 0 or 1 at the i'th entry.
%% Best values for parameters
tauF=6;%minuttes
Gmin=130.0;%mg/dL
Gmin3=1.5;%mg/dL/min
Gmin2=1.6;%mg/dL/min
%% Inisializing
k=length(Gm);

dGF=zeros(1,k);
%% preprocessing section
GFNS=NoiseSpikeFilter(Gm,dG);

%% Low pass filter

GF=lowpassfilter(GFNS,dt,tauF,k);

%% Lagrangian interpolation

dGF(1)=GF(1);
dGF(2)=GF(2);
% dGF(1)=0;
% dGF(2)=0;
for i=3:k
    T1=(t(i)-t(i-1))/...
        ((t(i-2)-t(i-1))*(t(i-2)-t(i)));
    T2=(t(i)-t(i-2))/...
        ((t(i-1)-t(i-2))*(t(i-1)-t(i)));
    T3=(2*t(i)-t(i-2)-t(i-1))/...
        ((t(i)-t(i-1))*(t(i)-t(i-2)));
    dGF(i)=T1*GF(i-2)...
        +T2*GF(i-1)...
        +T3*GF(i);
end

%% GRID
GRID=zeros(1,k);
for i=1:k
    %if GF(i)>Gmin && (max((dGF(end-2:end)>Gmin3)) || max((dGF(end-1:end)>Gmin2)))
    if GF(i)>Gmin && (max((dGF(i-2:i)>Gmin3))...
            || max((dGF(i-1:i)>Gmin2))) 
        GRID(i)=1;
    else
        %Not neccesary since GRID is a vector of zeros.
        GRID(i)=0;
    end
end



end