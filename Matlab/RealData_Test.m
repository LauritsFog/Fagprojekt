
%% Init
loadLib();
addpath('/Users/frederiknagel/Desktop/RealCGMData_DO_NOT_SHARE')

Days = 7;
InitData();

filename = 'Control-IQ_sample_Diasend.xls';
X = readtable(filename);

Gsc = table2array(X(:,2))*mmolL2mgdL;
Gsc = Gsc(1:Days*288);

%% Plot of data

T = 1:Days*288;
t = T*min2h;

figure(4);
plot(t,Gsc)
xlim([0 t(end)]);


%% Grid algo på real data

dg=10;
dt=10;
[GF,dGF,GRID_snack]=GridAlgo(Gsc,dg,dt,12,t);
x_snack=GRID_Filter(GRID_snack);

figure(4);
plot(t, Gsc,'b-',t,x_snack*350,'r-')
xlim([0 t(end)]);
ylabel({'CGM measurements', '[mg/dL]'});
xlabel('Time [h]');
title('Real Data - GRID alfo applied')





