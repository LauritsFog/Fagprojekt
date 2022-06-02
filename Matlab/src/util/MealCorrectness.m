function MealCorrectness(D,x,t)
%MEALCORRECTNESS is a function that computes how well the GRID algorithm
%has found meals based on the glucose concetration.

%   D: The mealplan with the correct times a meal is consumed
%   x: The output of the "GRID_Filter.m" algortihm. Its a vector with zeroes and ones
%      that indicates when the GRID algorithm finds a meal
%   t: The number of hours we allow the GRID algortihm to take to find a
%   meal


% --------- Basic info about mealplan and predicted meals ------------

a = nnz(D); % Finds number of meals in the mealplan
b = nnz(x); % finds number of predicted meals

X = sprintf('Number of meals:         %d',a);
Y = sprintf('Number of found meals:   %d',b);
disp(X)
disp(Y)


count = 0;

steps = t*12;

for i=1:length(x)
    if(D(i) > 0)    % tjek når der er et måltid
        
        for j=0:steps  % gå 1 time frem
            
            if(x(i+j)>0) % tjek for hvert step om der er et predicted måltid
               count = count+1;  % Hvis der er opdater
            end      
        end
        
    end    
end

P = round(count/b*100,2);
fprintf('Meals found within %d hour: %d \n' ,t,count);
X = sprintf('Procent:  %g\n',P);
disp(X)


counts = zeros(1,b);

id = 0;
for i=1:length(x)
    if(D(i) > 0) % tjek når der er et måltid
        id = id+1;
        j = 0; 
        while(x(i+j) == 0)
        j = j+1; 
        
        if i+j==length(D)
            break
        end
        
        end
        counts(1,id) = j;
        
    end    
end

Av = round((sum(counts)/b)*5,2);

X = sprintf('Average time to detect meal: %g min\n',Av);
disp(X)




end

