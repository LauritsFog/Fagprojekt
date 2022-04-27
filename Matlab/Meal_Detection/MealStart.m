function X=MealStart(G,t)
%% NOT FINISHED 
%%Finding the first time where a meal is detected

X=zeros(1,1);
for i=11:(length(G)-5)
if G(i)==1&&G(i-10)==0&&G(i-9)==0&&G(i-8)==0&&G(i-7)==0&&G(i-6)==0&&G(i-5)==0 ...
            &&G(i-4)==0&&G(i-3)==0&&G(i-2)==0&&G(i-1)==0&&G(i+1)==1&&G(i+2)==1&&G(i+3)==1 ...
            &&G(i+4)==1&&G(i+5)==1
        disp("YAY")
        X(end+1)=t(i);
end
end
        