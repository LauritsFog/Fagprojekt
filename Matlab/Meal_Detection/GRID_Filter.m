function [xTP, xFP]=GRID_Filter(D,GRID,t)
%%%                     Function to handle the Grid output
%Currently the GRID Algo returns a noisy output eg
%0,0,0,0,0,1,1,1,1,0,0,0,0,1,1,1,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1
%What GRID_Filter returns is vector of the same length but only with ones
%at the start of a meal eg.
%0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0
%% Input
%GRID the binary vector produced by GridAlgo
%% Output
%x the filtered binary vector
%% Merge meal peaks

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

%% The FIlter 
%The filter
% B=[0,1,1,1,1,1,1,1,1,1,1,1];
B=[0,1];
% Måske ikke det bedste filter, men det fanger sort set alle de store
% måltider
%
locationsOfMeal=strfind(GRID,B);
x=zeros(1,length(GRID));
% Adding 10 to get the position of the first 1. 
%(Just the number of 0 in the B vector before the first 1)
x([locationsOfMeal])=1;
% fjernet +10




%% Splits meals into True positive and False positive

xTP = x;
xFP = zeros(1,length(x));


for i=1:length(x)
   
    % Start if there is a predicted meal from the GRID_filter output
    if(x(i) == 1)
        T = 0;
        
        % For-loop for detecting if meal is a false postitive
        for j = 1:30 

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
            xTP(i) = 0;
            xFP(i) = 1;
            B = sprintf('False Positive at: %g',t(i));
            %disp(B)
        end
end

end

end
