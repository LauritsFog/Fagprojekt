function V = ComputeProcent(G)
%COMPUTEPROCENT Summary of this function goes here
%   This function computes how procent of the time the bloodsuger is in a
%   good value range??
%  - V = vector of the different procent
%  - The bloodsuger
%   

T = length(G);

% 54.0000   70.2000  180.0000  250.2000  664.0593

Great = 0; % the green part
Good = 0;% The Yellow part
bad = 0; % the orange part
critical = 0; % the red part




for i = 1:T
   
    %if bloodsugar is in the great part:
   if ge(G(i),70.2) && le(G(i),180.0)
       Great = Great+1;
   end
   
   %if bloodsugar is in the good part:
   if G(i)>180.0 && le(G(i),250.2)
       Good = Good+1;
   end
   
   % if bloodsugar is in the bad part:
   if G(i)<70.2 && G(i)>54.0
       bad = bad+1;
   end 
   if G(i)>250.0
       bad = bad+1;
   end 
   
   % if bloodsugar is in the critical part part:
   if le(G(i),54.0)
       critical = critical+1;
   end 
   
end


V = [(Great/T)*100 , (Good/T)*100, (bad/T)*100, (critical/T)*100];

                            % plot

% Initialize A                           
A = [];
A(1) = 100;
A(2) = 100-V(4);
A(3) = 100-V(4)-V(3);
A(4) = 100-V(4)-V(3)-V(2);


% Initialize the colors
Gcritcolors = {[255, 71, 71]/255; % Red
               [255, 154, 71]/255; % orange
               [255, 237, 71]/255; %Green
               [71, 255, 126]/255}; % yellow

       
figure(1)
for i = 1:1:length(A)
    if A(i)==0
        break;
    else
    area([0, 20],[A(i),A(i)],'FaceColor',Gcritcolors{i},'LineStyle','none')
    hold on
    end
end
plot(G);
xlim([0 20]);
ylim([0, 100]);
text(18,(A(1)-A(2))/2+A(2),"" + round(V(4),2) + "%")
text(18,(A(2)-A(3))/2+A(3),"" + round(V(3),2) + "%")
text(18,(A(3)-A(4))/2+A(4),"" + round(V(2),2) + "%")
text(18,(A(4))/2,"" + round(V(1),2) + "%")




end

