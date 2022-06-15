
%                 Test of GRID Algo on real data



%% Init
loadLib();
%addpath('/Users/frederiknagel/Desktop/RealCGMData_DO_NOT_SHARE')
addpath('G:\My Drive\Dtu\4 Semester\Fagprojekt')
Days = 3;
InitData();

filename = 'Control-IQ_sample_Diasend.xls';
X = readtable(filename);

Gsc = table2array(X(:,2))*mmolL2mgdL;
Gsc = Gsc(1:Days*288);

t = (1:1:Days*288)*5*min2h;


%% Meal vector

D = zeros(1,length(Gsc));
D(95) = 48;
D(188) = 45;
D(240) = 45;
D(256) = 45;
D(258) = 40;
D(287) = 55;

D(427) = 45;
D(470) = 45;
D(486) = 55;
D(521) = 55;
D(565) = 45;

D(699) = 45;
D(700) = 55;
D(746) = 50;
D(797) = 45;
D(848) = 45;

correctMeal = zeros(1,length(D));
for i=1:length(D)
   if D(i)>0
    correctMeal(i) = 1;
   end
end


%% Grid algo på real data

dg=3;
dt=1;
[GF,dGF,GRID]=GridAlgo(Gsc,dg,dt,[],t);

figure(4);
subplot(211)
for i = length(Gcrit):-1:1
    area([t0, tf]*min2h,[Gcrit(i),Gcrit(i)],'FaceColor',Gcritcolors{i},'LineStyle','none')
    hold on
end
plot(t, Gsc,'k-');
xlim([0 t(end)]);
ylim([0 350]);
ylabel({'CGM', '[mg/dL]'});
xlabel('Time [h]');
%title('Real Data - GRID algo')

subplot(212)
plot(t, Gsc,'k-',t,GRID*350,'r-',t,correctMeal*350,'g-','LineWidth',1.5)
xlim([0 t(end)]);
ylim([0 350]);
ylabel({'CGM', '[mg/dL]'});
xlabel('Time [h]');
legend({'CGM','Predicted Meal','Actual Meal'},'Position',[0.82 0.50 0.01 0.005])

%%
saveas(figure(4),[pwd '/Images/RealDataGRID.png']);

%%

tiledlayout(2,1)
% First plot
ax1 = nexttile;
plot(tspan,Gsc,'k',tspan,dGF,'b')
ylabel('CGM [mg/dL]')

% Second plot
ax2 = nexttile;
plot(tspan,dGF,'b',tspan,correctMeal*150,'g-')
title("Approximation G_F'")
xlabel('Time [min]')
ylabel("CGM' [mg/dl min]")




%%
%{
addpath('G:\My Drive\Dtu\4 Semester\Fagprojekt')
[testmealNoice, mealestNoice,dGF,ddGF]=MealSize(Gsc,t);

figure(2)
hold on
yyaxis right
plot(Gsc)
yyaxis left
plot(mealestNoice)
hold off



% Debug plot
tiledlayout(3,1)
% First plot
ax1 = nexttile;
yyaxis right
plot(t,Gsc,'k')
 yyaxis left
plot(t,mealestNoice)
title('Measurement noice')
xlabel('Time [min]')
ylabel('CGM [mg/dL]')

% Second plot
ax2 = nexttile;
plot(t,dGF,'b',t,correctMeal*150,'g-')
title("Approximation G_F'")
xlabel('Time [min]')
ylabel("CGM' [mg/dl min]")

% Third plot
ax3 = nexttile;
plot(t,ddGF,'g')
title("Approximation G_F''")
xlabel('Time [min]')
ylabel("CGM'' [mg/dl min^2]")
linkaxes([ax1 ax2 ax3],'x')

%}







