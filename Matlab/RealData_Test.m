
%                 Test of GRID Algo on real data



%% Init
loadLib();
addpath('/Users/frederiknagel/Desktop/RealCGMData_DO_NOT_SHARE')

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
[GF,dGF,GRID]=GridAlgo(Gsc,dg,dt,12,t);

figure(4);
hold on
plot(t, Gsc,'k-');
plot(t,GRID*350,'r-',t,correctMeal*350,'g-','LineWidth',1.5)
hold off
xlim([0 t(end)]);
ylim([0 350]);
ylabel({'CGM measurements', '[mg/dL]'});
xlabel('Time [h]');
legend('CGM','Predicted Meal','Actual Meal')
%title('Real Data - GRID algo')

saveas(figure(4),[pwd '/Images/RealDataGRID.png']);







