function [I] = get_image_name_OAM_1(Imagepath,current_time,channel,type)


I=imread(fullfile( Imagepath , ([ sprintf('img_000000%03d_', current_time) channel type])));



end


% nargout(@get_image_name_OAM_1)














