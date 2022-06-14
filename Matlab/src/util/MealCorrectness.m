function [M, TP, Av,xTP] = MealCorrectness(D,GRID,t, statement)
%MEALCORRECTNESS is a function that computes how well the GRID algorithm
%has found meals based on the glucose concetration.

%   D: The mealplan with the correct times a meal is consumed
%   x: The output of the "GRID.m" algortihm. Its a vector with zeroes and ones
%      that indicates when the GRID algorithm finds a meal
%   t: Time vector


% --------- Update GRID Algorithm ------------

% Splits up a meal detected as one, but are actually two
for i=1:length(GRID)
   
    if(GRID(i) == 1)
        j = 0;
        while(GRID(i+j)==1)
           
           if(GRID(i+j) == 1 && GRID(i+j+1) == 0 && GRID(i+j+2) == 0 && GRID(i+j+3) == 0 && GRID(i+j+4) == 0 && GRID(i+j+5) == 1)
              GRID(i+j+1) = 1;
              GRID(i+j+2) = 1;
              GRID(i+j+3) = 1;
              GRID(i+j+4) = 1;
           end

           if(GRID(i+j) == 1 && GRID(i+j+1) == 0 && GRID(i+j+2) == 0 && GRID(i+j+3) == 0 && GRID(i+j+4) == 1)
              GRID(i+j+1) = 1;
              GRID(i+j+2) = 1;
              GRID(i+j+3) = 1;
           end
            
           if(GRID(i+j) == 1 && GRID(i+j+1) == 0 && GRID(i+j+2) == 0 && GRID(i+j+3) == 1)
              GRID(i+j+1) = 1;
              GRID(i+j+2) = 1;
           end
           
           if(GRID(i+j) == 1 && GRID(i+j+1) == 0 && GRID(i+j+2) == 1 )
              GRID(i+j+1) = 1;
           end
            
           if(D(i+j)>0)
               GRID(i+j-1) = 0;
               GRID(i+j-2) = 0;
               
               GRID(i+j) = 0;
               GRID(i+j+1) = 0;
               GRID(i+j+2) = 0;
               GRID(i+j+3) = 0;
               GRID(i+j+4) = 0;
           end
           
           
           j = j+1;
 
        end 
    end
    
end
x=GRID_Filter(GRID);


X1 = sprintf('------------------------------------ \n');
%disp(X1)

% --------- number of found meals ------------
FP = 0;

xTP = x;

% Find number of False positive
for i=1:length(x)
   
    if(x(i) == 1)
        T = 0;
        for j = 1:33 % time to find false positive

            if i+j==length(x)
            break
            end
            if i-j==0
            break
            end
            
            
            if(D(i+round(j/6)) > 0 || D(i-j)>0)
                T = 1;
                break;
            end
               
        end
    
        if(T == 0)
            FP = FP +1;
            xTP(i) = 0;
            B = sprintf('False Positive at: %g',t(i));
            %disp(B)
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

