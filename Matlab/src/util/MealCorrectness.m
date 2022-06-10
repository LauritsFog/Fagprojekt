function [M, N, Av] = MealCorrectness(D,x)
%MEALCORRECTNESS is a function that computes how well the GRID algorithm
%has found meals based on the glucose concetration.

%   D: The mealplan with the correct times a meal is consumed
%   x: The output of the "GRID_Filter.m" algortihm. Its a vector with zeroes and ones
%      that indicates when the GRID algorithm finds a meal
%   t: The number of hours we allow the GRID algortihm to take to find a
%   meal


% --------- Basic info about mealplan and predicted meals ------------

M = nnz(D); % Finds number of meals in the mealplan
N = nnz(x); % finds number of predicted meals


% --------- Average time it take to find meal  ------------
counts = zeros(1,N);
id = 0;
for i=1:length(x)
    if(D(i) > 0) %  checks if there is a meal
        id = id+1;
        j = 0; 
        while(x(i+j) == 0)
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


%X = sprintf('Average time to detect meal: %g min\n',Av);
%disp(X)




end

