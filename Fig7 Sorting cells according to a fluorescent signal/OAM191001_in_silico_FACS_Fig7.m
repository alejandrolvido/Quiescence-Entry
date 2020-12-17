

%% Clustering code to distinguish mTFP1 positive from mTFP1 negative cells
%%%  In this example cells expressing CUP1p-Xbp1-mTFP1 and cells without it, 
%%%  undewent quiescence in the same microfluidic device, this script
%%%  sorts the cells into mTFP1(+) and mTFP1(-) for further analysis and
%%%  comparison.


% meaning of the columns in the Splot_info field:
% 1 cell number
% 2 position
% 3 cell exists
% 4 size at arrest
% 5 size at end
% 6 date
% 7 chamber
% 8 exp ID
% 9 Kmeans affiliation
% 10 Signal affiliation
% 11 ID on analysisY original
% 12 corrected chamber
% 13 corrected position

cd
clearvars
load('Toy_data_Fig6')

%% declare variables for plotting
        starti=1; % first time point
        endti=104; % last time point
        t1=starti:endti; % vector to define time window
        t=t1; % vector to define time fo plotting
        sr1=10; % sampling rate
        sr=12; % sampling rate conversion fro plotting in hours, 24 for every 5 minutes, 12 for every 10 min
        
%% declare matrix of the marker that will be used as reference for arrangements          
        lane_no=max(analysisY.Splot_info(:,7)); % number of lanes to be analyzed
        timeseriesMat=analysisY.SPLOT_6;% matrix containig the times series to be arranged        
       
%%%%%%%%%------- visualized initial data as Heat map ----------%%%%%%%%%%%%
      
        filteredmap1 = smoothActivityMap(timeseriesMat, 'SmoothParam', 0.9, 'UpSample', 1); % smooth function 
        figure(1);imagesc(filteredmap1);
        colorbar;colormap(jet)
        caxis([60 250])
        title('Initial Cdc10 Signal data')
        xlabel('Time(h)');ylabel('Cell Index')
        ax = gca;
        ax.XTickMode = 'manual';
        xticks(starti:sr:endti);
        curTick = ax.XTick;
        ax.XTickLabel = round((curTick)*sr1/60);
        
%%%%%%%%%%%%----------------  dendrogram -----------------%%%%%%%%%%%%%%%%%    
      
% The next steps will uses the Cdc10 matrix, which is contained in the
% structure "analysisY.SPLOT_6_F1" to arranged each time series as a
% dendrogram which is used for hierarchical clustering 

%%%%%%%%%%%%----------- loop for dendrogram --------------%%%%%%%%%%%%%%%%%    
       
        K = size(timeseriesMat, 1);
        Y = nan(K);
        for k = 1:K
            tmp = naneucdist(timeseriesMat(k, :), timeseriesMat);
            Y(:, k) = tmp; 
        end
        
 figure(2);
        tree2 = linkage(Y);  % option for hcl is 'single' (default)
        [~, ~, outperm2] = dendrogram(tree2, K);   
         title('Cdc10 Signal time series dendrogram')
        timeseriesMat_Cdc10= timeseriesMat(outperm2,t1);

% Arranged all markers according to the Cdc10 arrangement 'outperm2' 
analysisY.SPLOT_2_Cdc10=analysisY.SPLOT_2(outperm2,:);
analysisY.SPLOT_6_Cdc10=analysisY.SPLOT_6(outperm2,:);
analysisY.Splot_info_Cdc10=analysisY.Splot_info(outperm2,:);   

%% Inspect data as heat map to decide nunmber of clusters if required

%%%%%%----- visualize data according to dendrogram as heat map -----%%%%%%
        
        filteredmap2 = smoothActivityMap(timeseriesMat_Cdc10, 'SmoothParam', 0.9, 'UpSample', 1);
        figure(3); imagesc(filteredmap2);
        colorbar;colormap(jet)
        caxis([50 300])
        title(' Hierarchically arranged Cdc10 Signal time series')
        xlabel('Time(h)');ylabel('Cell Index')
        ax = gca;
        ax.XTickMode = 'manual';
        xticks(starti:sr:endti);
        curTick = ax.XTick;
        ax.XTickLabel = round((curTick)*sr1/60);


%%
%%%%%%%%%%--------------- k means clsutering ----------------%%%%%%%%%%%%%%    
    
% k-means, number of clustes (k) can be changed depending on inspection of the HCL 
% results, or other published methods to determine a relevant number of clusters (i.e. elbow method)
% however, once the clustering is validated using a second marker, the number of 
% clusters remain the same throughout the remaining experiments.  

%%%%%%%%%%--------------- k means clustering ----------------%%%%%%%%%%%%%%

Kgroup=3;% number of clusters
[grId,~,~] = kmeans(timeseriesMat_Cdc10, Kgroup, 'Replicates', 10, 'MaxIter', 10000, 'Display', 'final','Distance','correlation');
tabulate(grId) % 1st colum, cluster. 2nd column, number of cells, third column, percentage. 
analysisY.Splot_info_Cdc10(:,9)=grId; % asigns cluster affiliations to each time series in analysisY.Splot_info_Cdc10(:,9). 
groupedMat=cell(3,1);


%% ---- visualized clusters as a single heat map ----- %%

for ik=1:Kgroup % loop to concatenate the clusters on top of each other in a single heatmap
pooledMat_noNan = timeseriesMat_Cdc10((analysisY.Splot_info_Cdc10(:,9) == ik),:);
groupedMat{ik}= pooledMat_noNan;
end
 
groupedM=[groupedMat{1};groupedMat{2};groupedMat{3}]; % concatenation
filteredmap3 = smoothActivityMap(groupedM, 'SmoothParam', 0.9, 'UpSample', 1);
figure(4); imagesc(filteredmap3);
colorbar;colormap(jet)
caxis([50 300])
title('mixed mTFP1(+) and mTFP1(-) cells arranged into Cdc10-clusters')
xlabel('Time (h)');ylabel('Cell Index')
ax = gca;
ax.XTickMode = 'manual';
xticks(starti:sr:endti);
curTick = ax.XTick;
ax.XTickLabel = round((curTick)*sr1/60);

%% sorting of the control and experimental strains according to the presence of a mTFP1 signal  

% Transform potential irrationals into Nan before log transform
timeseriesM=analysisY.SPLOT_2_Cdc10; 
negNUM=find(timeseriesM<=0);
analysisY.timeseriesM_pos=timeseriesM;
analysisY.timeseriesM_pos(negNUM)=nan;
analysisY.timeseriesM_log=log10(analysisY.timeseriesM_pos);

% visualized log data as heat map
inputmap1 = analysisY.timeseriesM_log;
filteredmap1 = smoothActivityMap(inputmap1, 'SmoothParam', 0.9, 'UpSample', 1);
figure(5);imagesc(filteredmap1);
title('Initial log transformation of mTFP1 intensity')
colorbar;colormap(jet)

% imputation of NAN in each time series, imputed values assumed the average
% of the inmediately adjacent timepoints
noNAN_matrix=analysisY.timeseriesM_log;

for i=1:length(noNAN_matrix(:,1))
OUT=find(isnan(noNAN_matrix(i,:)));
if OUT~=0
for j=OUT
   if j==1
noNAN_matrix(i,j)=noNAN_matrix(i,j+1);
   elseif j==max(length(noNAN_matrix(1,:)))
noNAN_matrix(i,j)=noNAN_matrix(i,j-1);
   else
noNAN_matrix(i,j)=(noNAN_matrix(i,j-1)+noNAN_matrix(i,j+1))/2;
   end
end
else
end
end

% visualized imputed data
analysisY.timeseriesM_log_noNAN=noNAN_matrix;
inputmap2 = analysisY.timeseriesM_log_noNAN;
filteredmap2 = smoothActivityMap(inputmap2, 'SmoothParam', 0.9, 'UpSample', 1);
figure(6);imagesc(filteredmap2);
title('log(mTFP1 intensity), imputed')
colorbar;colormap(jet)

% collapse each time series to a value
AAA=1; %initial time point
BBB=104; % final time point

meanII=zeros(length(noNAN_matrix(:,1)),1);
for i= 1:length(noNAN_matrix(:,1))
meanI=nanmean(noNAN_matrix(i,AAA:BBB));
meanII(i,1)=meanI;
end

%%%%%% Kmeans on the meanI values for each time series, notice distance
%%%%%% metric is cityblock here. 

pooledMat_noNaN = meanII;
[grId,C,sumd] = kmeans(pooledMat_noNaN, 2, 'Replicates', 10, 'MaxIter', 10000, 'Display', 'final','Distance','cityblock');
tabulate(grId)

% asign cluster affiliation to ordered matrix;

analysisY.Splot_info_Cdc10(:,10)=grId;
G1=find(analysisY.Splot_info_Cdc10(:,10)==1);
G2=find(analysisY.Splot_info_Cdc10(:,10)==2);

% Asign cells with signal=1, cells without signal=2

G_1=nanmean(nanmean(meanII(G1)));
G_2=nanmean(nanmean(meanII(G2)));

    if G_1 <= G_2
analysisY.Splot_info_Cdc10(G2,10)=1;
analysisY.Splot_info_Cdc10(G1,10)=2;
    else  
    end

% sanity check to confirm that separation is correct.

cells_Signal = noNAN_matrix((analysisY.Splot_info_Cdc10(:,10)==1),:);
cells_No_signal = noNAN_matrix((analysisY.Splot_info_Cdc10(:,10)==2),:);
groupedMat = [cells_Signal;cells_No_signal];
inputmap3 = groupedMat;
filteredmap = smoothActivityMap(inputmap3, 'SmoothParam', 0.9, 'UpSample', 1);

figure(7); imagesc(filteredmap);
title('Sorted mTFP1(+) and mTFP1(-) cells')
colorbar;colormap(jet)
caxis([1 2.5])
figtmp.AlphaData = 1-isnan(inputmap3);
xlabel('Time (h)');ylabel('Cell index')
ax = gca;

ax.XTickMode = 'manual';
xticks(starti:sr:endti);
curTick = ax.XTick;
ax.XTickLabel = round((curTick)*sr1/60);
 
iso_cells=(analysisY.Splot_info_Cdc10(:,10)==2);
Fluo_cells=(analysisY.Splot_info_Cdc10(:,10)==1);

strain=1; % 1 to analysed the strain with the fluorescent signal, 2 to analyse the control strain 

if strain==1
    % declare matrices with cells having a Fluorescent signal
analysisY.SPLOT_2P=analysisY.SPLOT_2_Cdc10(Fluo_cells,:); 
analysisY.SPLOT_6P=analysisY.SPLOT_6_Cdc10(Fluo_cells,:); 
analysisY.Splot_infoP=analysisY.Splot_info_Cdc10(Fluo_cells,:);   
else
    
      % declare matrices with cells having not a Fluorescent signal
analysisY.SPLOT_2P=analysisY.SPLOT_2_Cdc10(iso_cells,:); 
analysisY.SPLOT_6P=analysisY.SPLOT_6_Cdc10(iso_cells,:); 
analysisY.Splot_infoP=analysisY.Splot_info_Cdc10(iso_cells,:); 
end

%%%%%%--------- obtain heatmaps and average curves for each marker------%%%

Kgroup=3; % define in previous step
colrs=['r'; 'b'; 'k'; 'g'; 'c';'m';'y']; % for plotting
marker=6; % marker to plot, 6=Cdc10.
limX=([0 endti]); 
group=(1:3);% specify number of clusters to plot
axis=[50 300]; % adjust according to marker intensity 

for ii=marker % loop through the channels to plot
      
%%% asign marker to analyse as Cdc10-arranged matrix
    
            if ii==1
                 timeseriesM=analysisY.SPLOT_1P;
                 titel1='Sfp1, Cdc10-clusters average curve';
            elseif ii==2
                 timeseriesM=analysisY.SPLOT_2P;
                 titel1='Stb3, Cdc10-clusters average curve';
            elseif ii==3
                 timeseriesM=analysisY.SPLOT_3P;
                 titel1='Xbp1, Cdc10-clusters average curve';
            elseif ii==4
                 timeseriesM=analysisY.SPLOT_4P;
                 titel1='Rtg1, Cdc10-clusters average curve';
            elseif ii==5
                 timeseriesM=analysisY.SPLOT_5P;
                 titel1='Gln3, Cdc10-clusters average curve';
            elseif ii==6
                 timeseriesM=analysisY.SPLOT_6P;
                 titel1='Cdc10-clusters average curve';
            end
            
            
timeseriesM(timeseriesM==0)=nan; % correction for cells that are born after time point 1
%%            
% calculate mean curves using each lane of the microfluidic device as
% biological replicate
numCells = 1:length(analysisY.Splot_infoP(:,7));
chamId = analysisY.Splot_infoP(:,7); 
CH_no=max(analysisY.Splot_infoP(:,7));
meanArr_perCham = cell(6,1);

for i=1:Kgroup % loop to obtain the average of each cluster per lane 
    meanArr_perCham{i} = nan(6, length(timeseriesM(1,:)));
    for j=1:CH_no
        ind = ((analysisY.Splot_infoP(:,9) == i) & (chamId == j));
        tmp = timeseriesM(ind,:);
        tmp1 = nanmean(tmp, 1);
        meanArr_perCham{i}(j, :) = tmp1;
    end
end

meanOfChambers = nan(Kgroup, length(meanArr_perCham{1}));
for i=1:Kgroup % loop to obtain the mean of each lane
    meanOfChambers(i,:) = nanmean(meanArr_perCham{i}, 1);
end
t=1:length(meanArr_perCham{1});


%%
f5=figure;

for ik4=group  % loop to plot the average curve for each cluster
    sem = std(meanArr_perCham{ik4}, [], 1, 'omitnan')./sqrt(size(meanArr_perCham{ik4},1));
    hold on
    s1 = shadedErrorBarV2(t, meanOfChambers(ik4,:), 2*sem, 'lineprops', colrs(ik4));
    hold on
end
             
    title(titel1)
    xlim(limX)
    xlabel('Time (h)');
    ax = gca;
    ax.XTickMode = 'manual';
    xticks(starti:sr:endti);
    curTick = ax.XTick;
    ax.XTickLabel = round((curTick)*sr1/60);

            if ii==1
                 titelx=[titel1 '_' num2str(ik4)];
                   ylabel('Nuclear intensity');
                   saveas(f5,titelx)
                   saveas(f5,'Sfp1_Kmeans','tif')                   
            elseif ii==2
                   ylabel('Nuclear intensity');                   
                   saveas(f5,'Stb3_Kmeans')
                   saveas(f5,'Stb3_Kmeans','tif')                   
            elseif ii==3
                   ylabel('Nuclear intensity');                   
                   saveas(f5,'Xbp1_Kmeans')
                   saveas(f5,'Xbp1_Kmeans','tif')
            elseif ii==4
                   ylabel('Nuclear intensity');                   
                   saveas(f5,'Rtg1_Kmeans')
                   saveas(f5,'Rtg1_Kmeans','tif')
            elseif ii==5
                   ylabel('Nuclear intensity');                   
                   saveas(f5,'Gln3_Kmeans')
                   saveas(f5,'Gln3_Kmeans','tif')
            elseif ii==6
                   ylabel('Cdc10 index');                   
                   saveas(f5,'Cdc10_Kmeans')
                   saveas(f5,'Cdc10_Kmeans','tif')
            end
            
 %%%%%%%%%%%%----------- obtain heat map ------------------%%%%%%%%%%%%%%%%  
 
groupedMat1=cell(3,1);
for ik=1:Kgroup % loop to concatenate the clusters on top of each other in a single heatmap
pooledMat_noNan = timeseriesM((analysisY.Splot_infoP(:,9) == ik),:);
groupedMat1{ik}= pooledMat_noNan;
end
groupedM1=[groupedMat1{1};groupedMat1{2};groupedMat1{3}];      

filteredmap = smoothActivityMap(groupedM1, 'SmoothParam', 0.9, 'UpSample', 1);

f6 = figure;imagesc(filteredmap);
colorbar;colormap(jet)
xlabel('Time (h)');ylabel('Cell Index')
ax = gca;
ax.XTickMode = 'manual';
xticks(starti:sr:endti);
curTick = ax.XTick;caxis(axis)
ax.XTickLabel = round((curTick)*sr1/60);
title(titel1)

            if ii==1
                title('Sfp1 nuclear translocation in Cdc10-Clusters')
                   saveas(f6,'Sfp1_heat')
                   saveas(f6,'Sfp1_heat','tif')
            elseif ii==2
                title('Stb3 nuclear translocation in Cdc10-Clusters')
                   saveas(f6,'Stb3_heat')
                   saveas(f6,'Stb3_heat','tif')
            elseif ii==3
                title('Xbp1 nuclear translocation in Cdc10-Clusters')
                   saveas(f6,'Xbp1_heat')
                   saveas(f6,'Xbp1_heat','tif')
            elseif ii==4
                title('Rtg1 nuclear translocation in Cdc10-Clusters')
                   saveas(f6,'Rtg1_heat')
                   saveas(f6,'Rtg1_heat','tif')
            elseif ii==5
                title('Gln3 nuclear translocation in Cdc10-Clusters')
                   saveas(f6,'Gln3_heat')
                   saveas(f6,'Gln3_heat','tif')
            elseif ii==6
                title('Cdc10 Signal Clusters')
                   saveas(f6,'Cdc10_heat')
                   saveas(f6,'Cdc10_heat','tif')
            end          
           
end 
   

