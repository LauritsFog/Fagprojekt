function PlotProcent(V)
 % plot

% Initialize A                           
A = [];
A(1) = 100;
A(2) = 100-V(1);
A(3) = 100-V(1)-V(2);
A(4) = 100-V(1)-V(2)-V(3);
A(5) = 100-V(1)-V(2)-V(3)-V(4);

% Initialize the colors
Gcritcolors = getCritColors;
       
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
    text(10,(A(5))/2,"" + round(V(5),2) + "%", 'FontSize', 18, 'HorizontalAlignment', 'center')
end

end


