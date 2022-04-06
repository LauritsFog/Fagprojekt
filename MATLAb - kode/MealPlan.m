function D = MealPlan(Days)
%MEALPLAN Summary of this function goes here
%   This function gives a random mealplan for a person




%Breakfast 
Bm = [90 100 110 120 130 140 150]/5;
Bt = [78:1:108];

% Lunch
Lm = [70 80 90 100]/5;
Lt = [138:1:162];

% Dinner
Dm = [70 80 90 100]/5;
Dt = [210:1:228];

n = Days;

% Meals: 3 meals a day for a month
D1 = zeros(288,1);
D = [];

for i = 1:n
   
    D1(randsample(Bt,1)) = randsample(Bm,1);
    D1(randsample(Lt,1)) = randsample(Lm,1);
    D1(randsample(Dt,1)) = randsample(Dm,1);
    
    D = [D; D1];
    D1 = zeros(288,1);
    
end


end

