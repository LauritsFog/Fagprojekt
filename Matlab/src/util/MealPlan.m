function D = MealPlan(Days,snack)
%MEALPLAN Summary of this function goes here
%   This function gives a random mealplan for a person

%{
Days = Number of days the mealplan should span over (Integer)
snack = If Mealplan should include snacks (Boolean)

Timeintevals for the different meals
Breakfast:   6.30 - 9.00
Lunch:       11.30 - 13.30
Dinner:      17.30 - 19.00
(If snack is included)
Snack 1:     15.00 - 16.00
Snack 2:     20.00 - 22.00
%}


%Breakfast 
Bm = [90 100 110 120 130 140 150]/5; %Sizes of meals
Bt = [78:1:108]; % Time inteval

% Lunch
Lm = [70 80 90 100]/5; % Sizes of meals
Lt = [138:1:162];      % Time inteval

% Dinner
Dm = [70 80 90 100]/5; % Sizes of meals
Dt = [210:1:228];       % Time inteval

% Snack
Sm = [18, 19, 20, 21, 22]/5; % Sizes of snacks
St1 = [180:1:192]; % time inteval for snack 1
St2 = [240:1:264]; % time inteval for snack 2

% Initialize a vector with zeroes that is 1 day long 
% (24h*60min / 5min = 1440/5 = 288)
% This means there is 288 5-min-intevals in a day
D1 = zeros(288,1);

%Initialize the complete mealplan vector
D = []; 

% Create the mealplan:
    for i = 1:Days
   
     % A random time and mealsize is chosen for Breakfast, Lunch and Dinner
     D1(randsample(Bt,1)) = randsample(Bm,1);
     D1(randsample(Lt,1)) = randsample(Lm,1);
     D1(randsample(Dt,1)) = randsample(Dm,1);
    
        % If snack is included a random time and snack size is chosen
        if snack == 1 
            D1(randsample(St1,1)) = randsample(Sm,1);
            D1(randsample(St2,1)) = randsample(Sm,1); 
        end
    
      % Add the day to the mealplan
      D = [D; D1];
    
      %Initialize the 1-day long vector again
      D1 = zeros(288,1);
    
    end


end

