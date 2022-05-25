function [V] = ComputeProcent(G, Gcrit)
%COMPUTEPROCENT Summary of this function goes here
%   This function computes how procent of the time the bloodsuger is in a
%   good value range??
%  - V = vector of the different procent
%  - The bloodsuger
%   

T = length(G);

% 54.0000   70.2000  180.0000  250.2000  664.0593

Optimalp2 = 0; % the orange part
Optimalp1 = 0;% The yellow part
Optimal = 0; % the green part
Optimalm1 = 0; % the pink part
Optimalm2 = 0; % the red part

for i = 1:T  
   % if bloodsugar is in the orange part:
   if Gcrit(4) < G(i) && G(i) < Gcrit(5)
       Optimalp2 = Optimalp2+1;
   end 
   
   %if bloodsugar is in the yellow part:
   if Gcrit(3) < G(i) && G(i) < Gcrit(4)
       Optimalp1 = Optimalp1+1;
   end
   
   % if bloodsugar is in the green part:
   if Gcrit(2) < G(i) && G(i) < Gcrit(3)
       Optimal = Optimal+1;
   end
   
   % if bloodsugar is in the pink part:
   if Gcrit(1) < G(i) && G(i) < Gcrit(2)
       Optimalm1 = Optimalm1+1;
   end
   
   % if bloodsugar is in the red part:
   if G(i) < Gcrit(1)
       Optimalm2 = Optimalm2+1;
   end
end


V = [(Optimalp2/T)*100, (Optimalp1/T)*100,(Optimal/T)*100, (Optimalm1/T)*100, (Optimalm2/T)*100];

                            % plot

% Initialize A                           
A = [];
A(1) = 100;
A(2) = 100-V(1);
A(3) = 100-V(1)-V(2);
A(4) = 100-V(1)-V(2)-V(3);
A(5) = 100-V(1)-V(2)-V(3)-V(4);

% Initialize the colors
Gcritcolors = {[255, 105, 105]/255;
               [255, 156, 156]/255;
               [156, 255, 159]/255;
               [255, 247, 156]/255;
               [255, 219, 156]/255};
       
for i = 1:1:length(A)
    if A(i)==0
        break;
    else
    area([0, 20],[A(i),A(i)],'FaceColor',Gcritcolors{length(A)-i+1},'LineStyle','none')
    hold on
    end
end
xlim([0 20]);
ylim([0, 100]);

for i = 1:length(A)-1
    if V(i) > 2
        text(10,(A(i)-A(i+1))/2+A(i+1),"" + round(V(i),2) + "%", 'FontSize', 18, 'HorizontalAlignment', 'center')
    end
end

if V(5) > 2
    text(10,(A(5))/2,"" + round(V(5),2) + "%")
end

end

