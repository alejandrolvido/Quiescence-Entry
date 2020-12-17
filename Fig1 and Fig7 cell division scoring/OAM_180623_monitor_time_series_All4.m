function[cell_data]=OAM_180623_monitor_time_series_All4(initial_tp,no_cell_of_interest,numbM,to_plot_here,colors_here,all_obj,zoom_in,mess,...
    exit_cond,Imagepath,fl_used,fl_used_2)

fl_used;
fl_used_2;
no_tp   =numbM;%
cell_data=-ones(1,4);
zoom_in;

disp(mess)

ccell = (all_obj.cells(:,:,initial_tp)==no_cell_of_interest);
ccell2= (all_obj.cells(:,:,initial_tp));

current_time=initial_tp;
OAM_180531_plot_subfigure_unvollstandiges4(to_plot_here,colors_here,no_cell_of_interest,current_time,ccell,zoom_in,ccell2,Imagepath,fl_used,fl_used_2,numbM);
              
checkbutton=1;

while checkbutton~=2
    [x1,y1,button] = ginput(1);
    if ~isempty(button)
        if (button==29) %forward in time
            current_time= current_time+1;
            if current_time>no_tp
                disp('end of time series')
            else
                ccell= (all_obj.cells(:,:,current_time)==no_cell_of_interest);
                ccell2= (all_obj.cells(:,:,current_time));
                OAM_180531_plot_subfigure_unvollstandiges4(to_plot_here,colors_here,no_cell_of_interest,current_time,ccell,zoom_in,ccell2,Imagepath,fl_used,fl_used_2,numbM);
            end
            %---------------------first +10 step------------------------------
        elseif button==57 %forward in time
            current_time= current_time+20;% 10
            if current_time>no_tp
                current_time=no_tp;
                disp('end of time series')
            else
                ccell= (all_obj.cells(:,:,current_time)==no_cell_of_interest);
                ccell2= (all_obj.cells(:,:,current_time));
                OAM_180531_plot_subfigure_unvollstandiges4(to_plot_here,colors_here,no_cell_of_interest,current_time,ccell,zoom_in,ccell2,Imagepath,fl_used,fl_used_2,numbM);
           
            end
            %---------------------end first +10 step--------------------------
        elseif  button==30 
            checkbutton=2;
        elseif (button==28)
            current_time= current_time-5;% 1
            if current_time<1
                disp('beginning of time series reached')
            else
                ccell= (all_obj.cells(:,:,current_time)==no_cell_of_interest);
                ccell2= (all_obj.cells(:,:,current_time));
                OAM_180531_plot_subfigure_unvollstandiges4(to_plot_here,colors_here,no_cell_of_interest,current_time,ccell,zoom_in,ccell2,Imagepath,fl_used,fl_used_2,numbM);
            end
            %--------------------- -10 step ---------------------
        elseif button==55
            current_time= current_time-20; % 10
            if current_time<1
                current_time=1;
                disp('beginning of time series reached')
            else
                ccell= (all_obj.cells(:,:,current_time)==no_cell_of_interest);
                ccell2= (all_obj.cells(:,:,current_time));
               OAM_180531_plot_subfigure_unvollstandiges4(to_plot_here,colors_here,no_cell_of_interest,current_time,ccell,zoom_in,ccell2,Imagepath,fl_used,fl_used_2,numbM);
            end
            %---------------end -10 step-------------------------------------
            elseif sum(button==exit_cond)>0
            cell_data(1:4)=[round(y1) round(x1) button current_time];
            checkbutton=2;
        else
            disp('odd error')
        end
    end
end

close all