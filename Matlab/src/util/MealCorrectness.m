function [M, TP, Av] = MealCorrectness(D,GRID, t, statement)
%MEALCORRECTNESS is a function that computes how well the GRID algorithm
%has found meals based on the glucose concetration.

%   D: The mealplan with the correct times a meal is consumed
%   x: The output of the "GRID.m" algortihm. Its a vector with zeroes and ones
%      that indicates when the GRID algorithm finds a meal
%   t: Time vector


% --------- Update GRID Algorithm ------------

% Splits up a meal detected as one, but are actually two


[xTP, xFP] = GRID_Filter(D,GRID,t);


X1 = sprintf('------------------------------------ \n');
%disp(X1)

% --------- number of found meals ------------

N = nnz(xTP)+nnz(xFP); % finds number of predicted meals
TP = nnz(xTP);
FP = nnz(xFP);

% --------- Basic info about mealplan and predicted meals ------------

M = nnz(D); % Finds number of meals in the mealplan

p = round(TP/M*100,2);

if(statement == 1)
X3 = sprintf('Actual Meals: %g',M);
X2 = sprintf('Total number of meals found: %g',N);
X1 = sprintf('Percent meals found: %g',p);
A = sprintf('Number of False Positive: %g',FP);
B = sprintf('Number of True Positive: %g',TP);
disp(X3)
disp(X2)
disp(A)
disp(B)
disp(X1)
end

% --------- Average time it take to find meal  ------------
counts = zeros(1,N);
id = 0;
for i=1:length(xTP)
    if(D(i) > 0) %  checks if there is a meal
        id = id+1;
        j = 0; 
        while(xTP(i+j) == 0)
        j = j+1;    % Antal 5-min tidsintervaller til næste predicted måltid
        if i+j==length(D)
            break
        end
        
        if j == 12*3
            j = 0;
            break
        end 
        
        end
        counts(1,id) = j;
        
    end    
end

Av = round((sum(counts)/nnz(counts))*5,2);

if(statement == 1)
X = sprintf('Average time to detect meal: %g min\n',Av);
 disp(X)
end



end

