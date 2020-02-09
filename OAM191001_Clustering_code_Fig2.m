
%% Clustering code to distinguish pattern of cell cycle arrest according to a Cdc10 marker


% Splot_info
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
load('Toy_data_Fig3')

% save('Toy_data_Fig3','analysisY')
        starti=1; % first time point
        endti=259; % last time point
        t1=starti:endti; % vector to define time window
        t=t1; % vector to define time fo plotting
        sr1=5; % sampling rate
        sr=24; % sampling rate conversion fro plotting in hours, 12 for every 5 minutes, 24 for every 10 min
        lane_no=max(analysisY.Splot_info(:,7)); % number of lanes to be analyzed
        timeseriesMat=analysisY.SPLOT_6;% matrix containig the times series to be arranged        
       
%%%%%%%%%------- visualized initial data as Heat map ----------%%%%%%%%%%%%
       
        filteredmap1 = smoothActivityMap(timeseriesMat, 'SmoothParam', 0.9, 'UpSample', 1); % smooth function 
        figure(1);imagesc(filteredmap1);
        colorbar;colormap(jet)
        caxis([60 250])
        title('Original Cdc10 Signal data')
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
 figure(2); title('dendrogram')
        tree2 = linkage(Y);  % option for hcl is 'single' (default)
        [~, ~, outperm2] = dendrogram(tree2, K);              
        timeseriesMat_Cdc10= timeseriesMat(outperm2,t1);
         
% Arranged all markers according to the Cdc10 arrangement 'outperm2' 
analysisY.SPLOT_6_Cdc10=analysisY.SPLOT_6(outperm2,:);
analysisY.Splot_info_Cdc10=analysisY.Splot_info(outperm2,:);   

%%
%%%%%%----- visualize data according to dendrogram as heat map -----%%%%%%
        
        filteredmap2 = smoothActivityMap(timeseriesMat_Cdc10, 'SmoothParam', 0.9, 'UpSample', 1);
        figure(3); imagesc(filteredmap2);
        colorbar;colormap(jet)
        caxis([50 300])
        title(' Arranged (HCL) Cdc10 Signal data')
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

%%%%%%%%%%--------------- clsutering using fcentroids----------------%%%%%%%%%%%%%%

Kgroup=3;% 
[grId,C,sumd] = kmeans(timeseriesMat_Cdc10, Kgroup, 'Replicates', 10, 'MaxIter', 10000, 'Display', 'final','Distance','correlation');
tabulate(grId) % 1st colum, cluster. 2nd column, number of cells, third column, percentage. 
analysisY.Splot_info_Cdc10(:,9)=grId; % asigns cluster affiliations to each time series in analysisY.Splot_info_Cdc10(:,9). 
groupedMat=cell(3,1);
for ik=1:Kgroup % loop to concatenate the clusters on top of each other in a single heatmap
pooledMat_noNan = timeseriesMat_Cdc10((analysisY.Splot_info_Cdc10(:,9) == ik),:);
groupedMat{ik}= pooledMat_noNan;
end
 
groupedM=[groupedMat{1};groupedMat{2};groupedMat{3}]; % correction  for concatenation

%%
%%%%%%----- visualized clusters as a single heat map -----%%%%%%

filteredmap3 = smoothActivityMap(groupedM, 'SmoothParam', 0.9, 'UpSample', 1);
figure(4); imagesc(filteredmap3);
colorbar;colormap(jet)
caxis([50 300])
title('k-means Cdc10-clusters')
xlabel('Time (h)');ylabel('Cell Index')
ax = gca;
ax.XTickMode = 'manual';
xticks(starti:sr:endti);
curTick = ax.XTick;
ax.XTickLabel = round((curTick)*sr1/60);

%%
%%%%%%%%%----- Get  mean curves and heat maps for each cluster ----%%%%%%%%

Kgroup=3; % define in previous step
colrs=['r'; 'b'; 'k'; 'g'; 'c';'m';'y']; % for plotting
marker=6; % marker to plot, chose between 1=Sfp1, 2=Stb3, 3=Xbp1, 4=Rtg1, 5=Gln3, 6=Cdc10.
limX=([0 endti]); 
group=(1:3);% specify number of clusters to plot
axis=[50 300]; % adjust according to marker intensity 

for ii=marker % loop through the channels to plot
      
%%% asign marker to analyse as Cdc10-arranged matrix
    
            if ii==1
                 timeseriesM=analysisY.SPLOT_1_Cdc10;
                 titel1='Sfp1, Cdc10-clusters average curve';
            elseif ii==2
                 timeseriesM=analysisY.SPLOT_2_Cdc10;
                 titel1='Stb3, Cdc10-clusters average curve';
            elseif ii==3
                 timeseriesM=analysisY.SPLOT_3_Cdc10;
                 titel1='Xbp1, Cdc10-clusters average curve';
            elseif ii==4
                 timeseriesM=analysisY.SPLOT_4_Cdc10;
                 titel1='Rtg1, Cdc10-clusters average curve';
            elseif ii==5
                 timeseriesM=analysisY.SPLOT_5_Cdc10;
                 titel1='Gln3, Cdc10-clusters average curve';
            elseif ii==6
                 timeseriesM=analysisY.SPLOT_6_Cdc10;
                 titel1='Cdc10-clusters average curve';
            end
            
%%            
% calculate mean curves using each lane of the microfluidic device as
% biological replicate
numCells = 1:length(analysisY.Splot_info_Cdc10(:,7));
chamId = analysisY.Splot_info_Cdc10(:,7); 
CH_no=max(analysisY.Splot_info_Cdc10(:,7));
meanArr_perCham = cell(6,1);

for i=1:Kgroup % loop to obtain the average of each cluster per lane 
    meanArr_perCham{i} = nan(6, length(timeseriesM(1,:)));
    for j=1:CH_no
        ind = ((analysisY.Splot_info_Cdc10(:,9) == i) & (chamId == j));
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
pooledMat_noNan = timeseriesM((analysisY.Splot_info_Cdc10(:,9) == ik),:);
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
   












