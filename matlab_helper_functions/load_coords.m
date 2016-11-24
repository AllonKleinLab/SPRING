function [xx,yy] = load_coords(coordinates_path)
	coords = csvread(coordinates_path);
    xx = coords(:,2);
    yy = coords(:,3);
end