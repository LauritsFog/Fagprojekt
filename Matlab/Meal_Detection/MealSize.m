function [M,MealEst]=MealSize(Gm,t,n)
%% Best values for parameters
tauF=6;%minuttes
dG=10;
dt=1;
%% Inisializing
%filtering to unique rows
% x=unique([t,Gm], 'rows');
% t=x(:,2);
% Gm=x(:,2);
%Making sure that Gm mod n is zero
%t=t(1:end-mod(length(Gm),n)+1);
%Gm=Gm(1:end-mod(length(Gm),n)+1);

%Extracting every n'th element from Gm
%Gm=Gm(1:n:end);
k=length(Gm);
dGF=zeros(1,k);
ddGF=zeros(1,k);
%% preprocessing section
GFNS=NoiseSpikeFilter(Gm,dG);
% Low pass filter
GF=lowpassfilter(GFNS,dt,tauF,k);

% Lagrangian interpolation

%dGF(1)=GF(1);
%dGF(2)=GF(2);
dGF(1)=0;
dGF(2)=0;
ddGF(1)=dGF(1);
ddGF(2)=dGF(2);
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
    %Second derivative
    ddGF(i)=T1*dGF(i-2)...
        +T2*dGF(i-1)...
        +T3*dGF(i);
end



%% Step 1
M=zeros(1,k);
for j=7:k
    for i =1:6
        if dGF(j-i)<0&& dGF(j-i+1)>0&& ddGF(j-i)<0&& ddGF(j-i+1)>0
            bool1=true;bool2=false;bool3=false;bool4=false;
            break;
        elseif i<=5&& dGF(j-i)<0.5&& dGF(j-i+1)>0.5&& ddGF(j-i)<0.005&& ddGF(j-i+1)>0.005
            bool1=false;bool2=true;bool3=false;bool4=false;
            break;
        elseif i<=5&& dGF(j-i)<1.25&& dGF(j-i+1)>1.25&& ddGF(j-i)<0.0125&& ddGF(j-i+1)>0.0125
            bool1=false;bool2=false;bool3=true;bool4=false;
            break;
        elseif i<=4&& dGF(j-i)<1.8&& dGF(j-i+1)>1.8&& ddGF(j-i)<0.018&& ddGF(j-i+1)>0.018
            bool1=false;bool2=false;bool3=false;bool4=true;
            break;
        else
            bool1=false;bool2=false;bool3=false;bool4=false;
        end
    end
    
    
    if bool1&& GF(j)>100&&dGF(j-1)<dGF(j)
        M(j)=1;
    elseif bool2&& GF(j)>100&&dGF(j-1)<dGF(j)
        M(j)=1.5;
    elseif bool3&& GF(j)>125&&dGF(j-1)<dGF(j)
        M(j)=2.25;
    elseif bool4&&GF(j)>120&&dGF(j-1)<dGF(j)
        M(j)=1.5;
    else
        M(j)=0;
    end
     
end

%% Step 2
for i=7:k
   
    if i-12<=0
        c=1;
    else
        c=i-12;
    end
    sumM=sum(M(c:i));
    
        greaterThanBool=sumM<=10;

    if greaterThanBool&&dGF(i-1)<ddGF(i-1)&& dGF(i)>ddGF(i)&& dGF(i-1)<dGF(i)&& ddGF(i-1)<ddGF(i)&& dGF(i)>0.5
        M(i)=3.0;
    elseif greaterThanBool&& ddGF(i-1)<dGF(i-1)&& ddGF(i)>dGF(i)&& dGF(i-1)<dGF(i)&& ddGF(i-1)<ddGF(i)&& dGF(i)>0.5
        M(i)=3.0;
    else
    end
        

end
%% Step 3
for i=7:k
   
    if i-6<=0
        c=1;
    else
        c=i-6;
    end
    sumM=sum(M(c:i));
    greaterThanBool=sumM<=12;

    if greaterThanBool&& M(i)==1.5 && GF(i)<300
        M(i)=3.5;
    elseif greaterThanBool&& M(i)==2.25 &&GF(i)<300
        M(i)=4.0;
    elseif 1<=sumM&&sumM<=5&&max(M(1:5))==1&&GF(i)<125
        M(i)=4.0;
    else
    end
        

end


%% Step 4
TDI=40;%U
for i=7:k
 

    if TDI<=28
        M(i)=0;
    elseif 45<=TDI&&TDI<=55&& GF(i)>130&&dGF(i)>0.15
        M(i)=1.75*M(i);
    elseif TDI>55&& GF(i)>130&& dGF(i)>0.15
        M(i)=3.25*M(i);
    else
    end
        

end


%% Step 5
for i=7:k
    if i-6<=0
        c=1;
    else
        c=i-6;
    end
    if sum(M(c:i))<=2
        M(c)=0;
    end
%     sumM=sum(M(c:i));
%     greaterThanBool=sumM<=15;
    
end


MealEst=4*M;
%% Extrastep













