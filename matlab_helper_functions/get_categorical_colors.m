function colors = get_categorical_colors(N)
% Convert a float between 0 and 1 into a hex color code using the jet
% colormap
% 
% ctracks  = Cell array. Each row is a color track. The first entry in the
%            row is the name of the track and subsequent entries give 
%            specific color values. If there are N cells, then this cell
%            array should have N+1 columns. 
% fname    = Path to save the color tracks file. 
%
%%
    colors = cell(N,1);
    cm = jet; 
    close(figure(1));
    for i=1:N
        frac = i/N;
        colorID = max(1, sum(frac > [0:1/length(cm(:,1)):1])); 
        rgb = cm(colorID, :);
        rgb = round(rgb*255);
        hex(:,2:7) = reshape(sprintf('%02X',rgb.'),6,[]).'; 
        hex(:,1) = '#';
        colors{i} = hex;
    end
end