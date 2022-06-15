function [M, TP, Av,xTP, FalsePositive] = MealCorrectness(D,GRID,t, statement)
%MEALCORRECTNESS is a function that computes how well the GRID algorithm
%has found meals based on the glucose concetration.

%   D: The mealplan with the correct times a meal is consumed
%   x: The output of the "GRID.m" algortihm. Its a vector with zeroes and ones
%      that indicates when the GRID algorithm finds a meal
%   t: Time vector


% --------- Update GRID Algorithm ------------

% Splits up a meal detected as one, but are actually two
for i=1:length(GRID)
   
    if(GRID(i) == 1) % Identifies when a meal is predicted
        j = 0;
        while(GRID(i+j)==1) % As long as the next index in the GRID oiutput is 1 this whileloop is going
           
            if(i+j+5 >= length(GRID)) % Makes sure we dont go out of bound
                break;
            end
            
            % If there is a 20min gap between two meals, its merged as one
            % meal
           if(GRID(i+j) == 1 && GRID(i+j+1) == 0 && GRID(i+j+2) == 0 && GRID(i+j+3) == 0 && GRID(i+j+4) == 0 && GRID(i+j+5) == 1)
              GRID(i+j+1) = 1;
              GRID(i+j+2) = 1;
              GRID(i+j+3) = 1;
              GRID(i+j+4) = 1;
           end

           % If there is a 15min gap between two meals, its merged as one
            % meal
           if(GRID(i+j) == 1 && GRID(i+j+1) == 0 && GRID(i+j+2) == 0 && GRID(i+j+3) == 0 && GRID(i+j+4) == 1)
              GRID(i+j+1) = 1;
              GRID(i+j+2) = 1;
              GRID(i+j+3) = 1;
           end
            
           % If there is a 10min gap between two meals, its merged as one
            % meal
           if(GRID(i+j) == 1 && GRID(i+j+1) == 0 && GRID(i+j+2) == 0 && GRID(i+j+3) == 1)
              GRID(i+j+1) = 1;
              GRID(i+j+2) = 1;
           end
           
           % If there is a 5min gap between two meals, its merged as one
            % meal
           if(GRID(i+j) == 1 && GRID(i+j+1) == 0 && GRID(i+j+2) == 1 )
              GRID(i+j+1) = 1;
           end
           
           % makes sure we dont go out of bound
           if(i+j-2 == 0)
               break;
           end
            
           % If there is an actual meal inside a row of 1's in the GRID
           % algorithm, it splits up the meal in two and makes some space
           
           if(D(i+j)>0)
               GRID(i+j-1) = 0;
               GRID(i+j-2) = 0;
               
               GRID(i+j) = 0;
               GRID(i+j+1) = 0;
               GRID(i+j+2) = 0;
               GRID(i+j+3) = 0;
               GRID(i+j+4) = 0;
           end
           
           % update value
           j = j+1;
 
        end 
    end
    
end
% apply GRID FIlter to the updated GRID output
x=GRID_Filter(GRID);


X1 = sprintf('------------------------------------ \n');
%disp(X1)

% --------- number of found meals ------------

%Initialize False Positive and Improved GRID filter
FP = 0;
xTP = x;
FalsePositive = zeros(1,length(x));

% Find number of False positive
for i=1:length(x)
   
    % Start if there is a predicted meal from the GRID_filter output
    if(x(i) == 1)
        T = 0;
        
        % For-loop for detecting if meal is a false postitive
        for j = 1:33 

            % Out of bound fix
            if i+j==length(x)
            break
            end
            if i-j==0
            break
            end
            
            % for each step checks if a real meal is within the timeinteval
            if(D(i+round(j/6)) > 0 || D(i-j)>0)
                T = 1;
                break;
            end
               
        end
    
        % If no meal is found its a false postitive
        if(T == 0)
            FP = FP +1;
            xTP(i) = 0;
            FalsePositive(i) = 1;
            B = sprintf('False Positive at: %g',t(i));
            disp(B)
        end
end

end


N = nnz(x); % finds number of predicted meals
TP = N - FP;

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

