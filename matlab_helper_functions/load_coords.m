function [xx,yy] = load_coords(coordinates_path)
% Load coordinates that have been saved by SPRING
%
% coordinates_path = Path to the coordinates file
% xx               = Array of x-coordinates
% yy               = Array of y-coordinates
%
%%
    coords = csvread(coordinates_path);
    xx = coords(:,2);
    yy = coords(:,3);
end