function hex = frac_to_hex(frac)
    cm = jet; 
    close(figure(1));
    colorID = max(1, sum(frac > [0:1/length(cm(:,1)):1])); 
    rgb = cm(colorID, :);
    rgb = round(rgb*255);
    hex(:,2:7) = reshape(sprintf('%02X',rgb.'),6,[]).'; 
    hex(:,1) = '#';
end