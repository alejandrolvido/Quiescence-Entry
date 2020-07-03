
%%% this script is an example of how to mark the time of cell division by
%%% looking at the Cdc10 signal (subplot top left) and budding (subplot bottom left)
%%% and other two fluorescent channels (subplots on the right)
pwd
currentFolder = pwd;

exp_name='OAM_180915_6c_two_q_entries_120ms';
pos_here=3; %
Path_1=[currentFolder '\'];% 'C:\Users\aleja\Documents\MATLAB\'; 

%%% loads images which have been unzipped in the folder 'Fig1 and Fig7 cell division scoring'
%%% if the unzipping process creates a Pos3 folder within the Pos3 folder
%%% use the next line otherwise used the next one.
Imagepath = [ Path_1 exp_name '\Pos' num2str(pos_here) '\Pos' num2str(pos_here) '\']; % if unzipped creates Pos3 folder within Pos3 folder 
% Imagepath = [ Path_1 exp_name '\Pos' num2str(pos_here) '\']; % if unzipped creates only one Pos3 folder


% % %loads segmentation
seg_file=[Path_1 exp_name '\' exp_name '_pos' num2str(pos_here) 'dm1'];
load(seg_file);

% %loads extracted fluorescent data
fl_data=[Path_1 exp_name '\' exp_name '_pos' num2str(pos_here) '_extr_new2'];
load(fl_data); 

% load scores with cells status
sta_file=[Path_1 exp_name '\' exp_name '_pos' num2str(pos_here) '_cell_status']; %windows
load(sta_file); 
   
%%% structure to store data
pos_data.cell_number        = zeros(1,no_obj); % yep, this is the number
pos_data.pos_number         = zeros(1,no_obj); % just the position
pos_data.Cdc10_Onset        = zeros(5,no_obj); % Spc42 separation time 

% select fluorescent channels
% select FL 1: TFP GFP/mNG mKok mRuby mNeptune CYOFP
fl_used=[0 0 0 0 0 1 0];
% select FL 2: TFP GFP/mNG mKok mRuby mNeptune CYOFP
fl_used_2=[1 0 0 0 0 0 0];

F1=1; % time point to start
F2=30; % end time point
Fcells1=find(cell_exists(:,2)<=F1); % cells that exits before the return to rich conditions
initial_tp=20;%%
numbM=F2;

for  i=Fcells1'%1:length(all_obj.Cdc10_STD(:,1))
   
    disp(['cell no = ' num2str(i) '/' num2str(no_obj)])
    % ------ assignments -------
    pos_data.cell_number(1,i)=i;
    pos_data.pos_number(1,i) =pos_here;

    W5=all_obj.Cdc10_STD(i,1:F2); % plot Cdc10 Signal
    
    to_plot_here=[W5];
    colors_here=[1 0.6 0;1 0 0;0 0 0; 0 0 1; 0 1 0]; %
    no_cell_of_interest=i;
    zoom_in=[1 40];%second parameter is cell margin
    
       if pos_status.cell_status(1,i)==1
              exit_cond=32;
        mess=(['cell no: ' num2str(no_cell_of_interest) '/' num2str(no_obj) ', mark first division']);
        [cell_data]=OAM_180623_monitor_time_series_All4(initial_tp,no_cell_of_interest,numbM,to_plot_here,colors_here,all_obj,zoom_in,mess,...
    exit_cond,Imagepath,fl_used,fl_used_2);
       pos_data.Cdc10_Onset(1,i)=cell_data(1,4); 
       else
       end
       
       if  pos_data.Cdc10_Onset(1,i)>0
       exit_cond=32;
              initial_tp=pos_data.Cdc10_Onset(1,i);
        mess=(['cell no: ' num2str(no_cell_of_interest) '/' num2str(no_obj) ',  mark second division']);
        [cell_data]=OAM_180623_monitor_time_series_All4(initial_tp,no_cell_of_interest,numbM,to_plot_here,colors_here,all_obj,zoom_in,mess,...
    exit_cond,Imagepath,fl_used,fl_used_2);
       pos_data.Cdc10_Onset(2,i)=cell_data(1,4); 
       else 
       end
       
        if  pos_data.Cdc10_Onset(2,i)>0
         exit_cond=32;
               initial_tp= pos_data.Cdc10_Onset(2,i);
        mess=(['cell no: ' num2str(no_cell_of_interest) '/' num2str(no_obj) ', mark third division']);
        [cell_data]=OAM_180623_monitor_time_series_All4(initial_tp,no_cell_of_interest,numbM,to_plot_here,colors_here,all_obj,zoom_in,mess,...
    exit_cond,Imagepath,fl_used,fl_used_2);
       pos_data.Cdc10_Onset(3,i)=cell_data(1,4); 
        else 
        end
       
         if  pos_data.Cdc10_Onset(3,i)>0
        exit_cond=32;
        initial_tp=pos_data.Cdc10_Onset(3,i);
        mess=(['cell no: ' num2str(no_cell_of_interest) '/' num2str(no_obj) ',  mark fourth division']);
        [cell_data]=OAM_180623_monitor_time_series_All4(initial_tp,no_cell_of_interest,numbM,to_plot_here,colors_here,all_obj,zoom_in,mess,...
    exit_cond,Imagepath,fl_used,fl_used_2);
       pos_data.Cdc10_Onset(4,i)=cell_data(1,4); 
        else 
         end
        
       if  pos_data.Cdc10_Onset(4,i)>0
       exit_cond=32;
        initial_tp=pos_data.Cdc10_Onset(4,i);
        mess=(['cell no: ' num2str(no_cell_of_interest) '/' num2str(no_obj) ',  mark fifth division']);
        [cell_data]=OAM_180623_monitor_time_series_All4(initial_tp,no_cell_of_interest,numbM,to_plot_here,colors_here,all_obj,zoom_in,mess,...
    exit_cond,Imagepath,fl_used,fl_used_2);
       pos_data.Cdc10_Onset(5,i)=cell_data(1,4); 
       else
       end
       
    disp('======================================================================')
    disp(['pos_data.cell_number     = ' num2str(pos_data.cell_number(1,i))])
    disp(['pos_data.pos_number      = ' num2str(pos_data.pos_number(1,i))])
    disp(['pos_data.Cdc10_Onset =' num2str(pos_data.Cdc10_Onset(1,i))])

         name1=[exp_name '_pos' num2str(pos_here) '_Exit_Cdc10'];
%          save(name1, 'pos_data','pos_here')

 end
 
 
