function [x_cn,y_cn]=get_wind_coord1(Ihere,cell_margin)

[x_size,y_size]=size(Ihere);

if sum(Ihere(:))>0 % not empty object
s_f=regionprops(Ihere,'Centroid','BoundingBox');
if ~isempty(s_f)   
    bbox=round(s_f(1).BoundingBox);
    lower_x_limit=max(1,bbox(1)-cell_margin);
    upper_x_limit=min(y_size,bbox(1)+bbox(3)+cell_margin);
    lower_y_limit=max(1,bbox(2)-cell_margin);
    upper_y_limit=min(x_size,bbox(2)+bbox(4)+cell_margin);
    x_cn=lower_x_limit:upper_x_limit;
    y_cn=lower_y_limit:upper_y_limit;
else
    disp('empty or multiple object given - error')
end

else
    x_cn=1:y_size;
    y_cn=1:x_size;
end