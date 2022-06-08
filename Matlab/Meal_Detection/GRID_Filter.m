function x=GRID_Filter(GRID)
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
%%
%The filter
B=[0,1,1,1,1,1,1,1,1];
% Måske ikke det bedste filter, men det fanger sort set alle de store
% måltider

%
locationsOfMeal=strfind(GRID,B);
x=zeros(1,length(GRID));
% Adding 10 to get the position of the first 1. 
%(Just the number of 0 in the B vector before the first 1)
x([locationsOfMeal])=1;
% fjernet +10
end
