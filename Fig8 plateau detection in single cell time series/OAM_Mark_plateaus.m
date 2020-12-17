
%%% this script loads the time series from experiments with two successive
%%% quiescence entries and calculates the time when the nuclear
%%% accumulation of stress transcription factors reaches a plateau
%%% The automatic asigments of plateau times are then manually curated to
%%% exclude invalid asigments. Finally the plateau times from the 1st Q are
%%% compare to the plateau times from the 2nd Q


%% load analysis type MasterM 
load('Stb3_time_series') 

MasterM=Stb3_time_series;
MasterM2=cell(size(MasterM));
MasterMB=cell(size(MasterM));
F1=70;
F2=350;
nmax=numel(MasterM2);
%cellsA=find(MasterM{nmax}(:,3)>=F1 & MasterM{nmax}(:,3)<=F2 ); 
cellsA=find(MasterM{nmax}(:,3)<=F1);
cellsB=find(MasterM{nmax}(:,3)<=F2);

% get cells present during 1st Q
for ik=1:nmax
MasterM2{ik}=MasterM{ik}(cellsA,:);
end


% get cells present during 2nd Q
for ik=1:nmax
MasterMB{ik}=MasterM{ik}(cellsB,:);
end

nfactor=8;% compression factor
Ns=1; % in this case 1 because Stb3 is the first channel
Mat1=MasterM2{Ns};
[Mat_compress,GR_Matrix] = find_maxgrowthrate(Mat1, nfactor);

AAA=round(70/nfactor);
BBB=round(240/nfactor);
plateau_t=zeros(3,length(GR_Matrix(:,1)));

for iy=1:length(GR_Matrix(:,1))
    % Find out where the standard deviation is less than 0.2
    index = find(GR_Matrix(iy,AAA:BBB)==max(GR_Matrix(iy,AAA:BBB)));
    index2 = find(GR_Matrix(iy,AAA+index:BBB)<0,1,'first');
        
    if numel(index)>=2 ||  isempty(index2)
      plateau_t(:,iy)=NaN; 
    else
      plateau_t(1,iy)=index;
      index2=index+index2;
      plateau_t(2,iy)=index2;
      plateau_t(3,iy)=(round((index+index2)/2));
    end
end

plateau_t=(plateau_t(3,:)+AAA)*nfactor;
nanmean(plateau_t)
nanstd(plateau_t)

% nfactor=16;
Mat1=MasterMB{Ns};
[Mat_compress,GR_Matrix] = find_maxgrowthrate(Mat1, nfactor);

AAA1=round(350/nfactor);
BBB1=round(483/nfactor)-1;
plateau_t2=zeros(1,length(GR_Matrix(:,1)));

for iy=1:length(GR_Matrix(:,1))
    % Find out where the standard deviation is less than 0.2
    index = find(GR_Matrix(iy,AAA1:BBB1)==max(GR_Matrix(iy,AAA1:BBB1)));
    index2 = find(GR_Matrix(iy,AAA1+index:BBB1)<0,1,'first');
        
    if numel(index)>=2 ||  isempty(index2)
      plateau_t2(:,iy)=NaN; 
    else
      index2=index+index2;
      plateau_t2(1,iy)=(round((index+index2)/2));
    end
end

plateau_t2=(plateau_t2(1,:)+AAA1)*nfactor;
nanmean(plateau_t2)
nanstd(plateau_t2)

%% quality check by inspection of the plateau asigments for each cell

% First entry quality check

%%% this loops show the time series with vertical dotted lines that mark 
%%% the changes of media and the asigments the for start of plateaus in the time series. 
%%% spacebar to accept the asigment as valid, mouse left click to reject it
%%% as invalid. red line=plateau, blue line=starvation, green line=return
%%% to growth

arrestT=70;
arrestT2=250;
returnT=350;
Mat_C=plateau_t;
cell_data=zeros(1,size(Mat_C,2));

 for iy=1:size(Mat_C,2)
     
   AA1=Mat_C(1,iy);
   if ~isnan(AA1)
    
   figure
     plot((MasterM2{Ns}(iy,:)));
      xline(AA1,'--r') ;%  'Starvation' 
      xline(arrestT,'--b') ;%  'Starvation' 
      xline(returnT,'--b') ;%  'Starvation' 
      xline(arrestT2,'--g') ;% 
      
      title(['cell no= ' num2str(iy) '  spacebar (accept) or left click (reject)']);
    
     exit_cond=32; %space 
     
     [x1,y1,button] = ginput(1);
    
        if button==exit_cond
            cell_data(1,iy)=1;
            close
        else
           close 
        end
 
   else
   end
 end
 
 okcells=find(cell_data==1);
 
 Plat1=(Mat_C(1,okcells)-arrestT);
 nanmean(Plat1)
 nanstd(Plat1)
 
 
 % Second entry quality check
 
%%% this loops show the time series with vertical dotted lines that mark 
%%% the changes of media and the asigments the for start of plateaus in the time series. 
%%% spacebar to accept the asigment as valid, mouse left click to reject it
%%% as invalid.

arrestT=70;
arrestT2=350;
returnT=250;
Mat_C=plateau_t2;
cell_data=zeros(1,size(Mat_C,2));

 for iy=1:size(Mat_C,2)
     
   AA1=Mat_C(1,iy);
   if ~isnan(AA1)
    
   figure
     plot((MasterMB{Ns}(iy,:)));
      xline(AA1,'--r') ;%  'Starvation' 
      xline(arrestT,'--b') ;%  'Starvation' 
      xline(returnT,'--b') ;%  'Starvation' 
      xline(arrestT2,'--g') ;% 
      
      title(['cell no= ' num2str(iy) '  spacebar (accept) or left click (reject)']);
    
     exit_cond=32; %space 
     
     [x1,y1,button] = ginput(1);
    
        if button==exit_cond
            %0 for mother 1 for daughter
            cell_data(1,iy)=1;
            close
        else
           close 
        end
 
   else
   end
 end
 
 okcells2=find(cell_data==1);
 
 Plat2=(Mat_C(1,okcells2)-arrestT2);
 nanmean(Plat2)
 nanstd(Plat2)
 
%% creates matrix with plateau times for 1st and 2nd Q entry after all plateau asigments have been verified 
 
Plats=nan(size(Plat2,2),2);
Plats(1:size(Plat1,2),1)=Plat1';
Plats(1:size(Plat2,2),2)=Plat2';

Platsh=Plats*5/60;

boxplot(Platsh)

title('Time to reach maximum nuclear accumulation');
xticklabels({'1st Q','2nd Q'});
xlabel('Stress Transcription factors');
ylabel('hours');

titel1='Stb3'; % name after the channel your looking at
path_save='C:\Users\aleja\OneDrive\Documents\MATLAB doc\';
export_path=[path_save titel1];
hgexport(gcf, export_path, hgexport('factorystyle'), 'Format', 'pdf');

%% statistical test for differences between 1st and 2nd Q plateau times.

 [K,P]=kstest2(Platsh(:,1),Platsh(:,2));

 save('Stb3_Plat','Platsh','K','P','Ns','arrestT','arrestT2','returnT'); 

