
%% Clustering code to distinguish cells that crossed the START checkpoint before Starvation onset


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
load('Toy_data_Fig4');

%% Cdc10 Arrangement 

        starti=1; % first time point
        endti=214; % last time point
        t1=starti:endti; % vector to define time window
        t=t1; % vector to define time fo plotting
        sr1=5; % sampling rate
        sr=24; % sampling rate conversion fro plotting in hours, 12 for every 5 minutes, 24 for every 10 min
        lane_no=max(analysisZ.Splot_info(:,7)); % number of lanes to be analyzed
        timeseriesMat=analysisZ.SPLOT_6;% matrix containig the times series to be arranged             
       
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
analysisZ.SPLOT_6_Cdc10=analysisZ.SPLOT_6(outperm2,:);
analysisZ.Splot_info_Cdc10=analysisZ.Splot_info(outperm2,:);

%%
%%%%%%----- visualize data according to dendrogram as heat map -----%%%%%%
        
        filteredmap2 = smoothActivityMap(timeseriesMat_Cdc10, 'SmoothParam', 0.9, 'UpSample', 1);
        figure(1); imagesc(filteredmap2);
        colorbar;colormap(jet)
        caxis([50 300])
        title(' Hierarchically Arranged Cdc10 Signal data')
        xlabel('Time(h)');ylabel('Cell Index')
        ax = gca;
        ax.XTickMode = 'manual';
        xticks(starti:sr:endti);
        curTick = ax.XTick;
        ax.XTickLabel = round((curTick)*sr1/60);


%% Kmeans for Cdc10 on the windown of time to define cells that crossed START ebfore starvation onset
AAA=12; % start of starvation on lane
BBB=80; % end of evaluated windown of time 
Mat_AR=analysisZ.SPLOT_6_Cdc10;
%% collapse time series during the windown of time to a mean value
meanII=zeros(length(analysisZ.SPLOT_6_Cdc10(:,1)),1);
for i= 1:length(analysisZ.SPLOT_6_Cdc10(:,1))
meanI=nanmean( Mat_AR(i,AAA:BBB));
meanII(i,1)=meanI;
end

%%%%%% Kmeans on specified windown of time, notice used of cityblock
%%%%%% as distance metric

pooledMat_noNaN =meanII;
[grId,~,~] = kmeans(pooledMat_noNaN, 2, 'Replicates', 10, 'MaxIter', 10000, 'Display', 'final','Distance','cityblock');
tabulate(grId)

%asign group to ordered matrix;
analysisZ.Splot_info_Cdc10(:,10)=grId;
G1=find(analysisZ.Splot_info_Cdc10(:,10)==1);
G2=find(analysisZ.Splot_info_Cdc10(:,10)==2);
G_1=nanmean(nanmean(meanII(G1)));
G_2=nanmean(nanmean(meanII(G2)));

    if G_1 <= G_2
analysisZ.Splot_info_Cdc10(G2,10)=1; % cell that crossed START shortly before starvation onset
analysisZ.Splot_info_Cdc10(G1,10)=2; % cell that crossed START at other time or never crossed START after starvation onset
    else  
    end

    % sanity check, visualization of separated populations, START-crossers on top  
pooledMat_noNan_g1 = Mat_AR((analysisZ.Splot_info_Cdc10(:,10)==1),:);
pooledMat_noNan_g2 = Mat_AR((analysisZ.Splot_info_Cdc10(:,10)==2),:);
groupedMat = [pooledMat_noNan_g1;pooledMat_noNan_g2];
inputmap = groupedMat;
filteredmap = smoothActivityMap(inputmap, 'SmoothParam', 0.9, 'UpSample', 1);
figure(2);imagesc(filteredmap);
colorbar;colormap(jet)
title('Sorting of cells according START crossing at starvation onset')
xlabel('Time (h)');ylabel('Cell Index')

START_NOT_CROSSED_BEFORE=find(analysisZ.Splot_info_Cdc10(:,10)==2)';
START_CROSSED_BEFORE=find(analysisZ.Splot_info_Cdc10(:,10)==1)';

strain=1; % 1 for START-crossers
if strain==1
% declare matrices with cells having crossed START at the onset of
% starvation
analysisZ.SPLOT_6P=analysisZ.SPLOT_6_Cdc10(START_CROSSED_BEFORE,:); 
analysisZ.Splot_infoP=analysisZ.Splot_info_Cdc10(START_CROSSED_BEFORE,:);   
else
% declare matrices with cells having crossed START at other times
analysisZ.SPLOT_6P=analysisZ.SPLOT_6_Cdc10(START_NOT_CROSSED_BEFORE,:); 
analysisZ.Splot_infoP=analysisZ.Splot_info_Cdc10(START_NOT_CROSSED_BEFORE,:); 
end


%%
%%%%%%%%%%--------------- k means clsutering ----------------%%%%%%%%%%%%%%    
    
% k-means, number of clustes (k) can be changed depending on inspection of the HCL 
% results, or other published methods to determine a relevant number of clusters (i.e. elbow method)
% however, once the clustering is validated using a second marker(budding in our case), the number of 
% clusters remain the same throughout the remaining experiments.  

%%%%%%%%%%--------------- k means clsutering ----------------%%%%%%%%%%%%%%
timeseriesMat2=analysisZ.SPLOT_6P; % matrix of marker used to establish clusters 
Kgroup=2;% number of cluster 
[grId,~,~] = kmeans(timeseriesMat2, Kgroup, 'Replicates', 10, 'MaxIter', 10000, 'Display', 'final','Distance','correlation');
tabulate(grId) % 1st colum, cluster. 2nd column, number of cells, third column, percentage. 
analysisZ.Splot_infoP(:,9)=grId; % asigns cluster affiliations to each time series in analysisZ.Splot_info_Cdc10(:,9). 
groupedMat=cell(5,1);
for ik=1:Kgroup % loop to concatenate the clusters on top of each other in a single heatmap
pooledMat_noNan = timeseriesMat2((analysisZ.Splot_infoP(:,9) == ik),:);
groupedMat{ik}= pooledMat_noNan;
end

groupedM=[groupedMat{1};groupedMat{2};groupedMat{3};groupedMat{4};groupedMat{5}]; % custom concatenation 

%%
%%%%%%----- visualized clusters as a single heat map -----%%%%%%

filteredmap3 = smoothActivityMap(groupedM, 'SmoothParam', 0.9, 'UpSample', 1);
figure(3); imagesc(filteredmap3);
colorbar;colormap(jet)
caxis([50 350])
title('Cdc10-clusters of cells that crossed START at starvation onset')
xlabel('Time (h)');ylabel('Cell Index')
ax = gca;
ax.XTickMode = 'manual';
xticks(starti:sr:endti);
curTick = ax.XTick;
ax.XTickLabel = round((curTick)*sr1/60);

%%%%%%--------- obtain heatmaps and average curves for each marker------%%%

Kgroup=2; % define in previous step
colrs=['r'; 'b'; 'k'; 'g'; 'c';'m';'y']; % for plotting
marker=6; % marker to plot, 6=Cdc10.
limX=([0 endti]); 
group=(1:2);% specify number of clusters to plot
axis=[50 350]; % adjust according to marker intensity 

for ii=marker % loop through the channels to plot
      
%%% asign marker to analyse as Cdc10-arranged matrix
    
            if ii==1
                 timeseriesM=analysisZ.SPLOT_1P;
                 titel1='Sfp1, Cdc10-clusters average curve';
            elseif ii==2
                 timeseriesM=analysisZ.SPLOT_2P;
                 titel1='Stb3, Cdc10-clusters average curve';
            elseif ii==3
                 timeseriesM=analysisZ.SPLOT_3P;
                 titel1='Xbp1, Cdc10-clusters average curve';
            elseif ii==4
                 timeseriesM=analysisZ.SPLOT_4P;
                 titel1='Rtg1, Cdc10-clusters average curve';
            elseif ii==5
                 timeseriesM=analysisZ.SPLOT_5P;
                 titel1='Gln3, Cdc10-clusters average curve';
            elseif ii==6
                 timeseriesM=analysisZ.SPLOT_6P;
                 titel1='Cdc10-clusters average curve';
            end
            
%%            
% calculate mean curves using each lane of the microfluidic device as
% biological replicate
numCells = 1:length(analysisZ.Splot_infoP(:,7));
chamId = analysisZ.Splot_infoP(:,7); 
CH_no=max(analysisZ.Splot_infoP(:,7));
meanArr_perCham = cell(6,1);

for i=1:Kgroup % loop to obtain the average of each cluster per lane 
    meanArr_perCham{i} = nan(6, length(timeseriesM(1,:)));
    for j=1:CH_no
        ind = ((analysisZ.Splot_infoP(:,9) == i) & (chamId == j));
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
             
    title('Non G1 (red) and Post mitotic G1 (blue) Q cells');
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
                    saveas(f5,'Sfp1_Kmeans','pdf')  
            elseif ii==2
                   ylabel('Nuclear intensity');                   
                   saveas(f5,'Stb3_Kmeans')
                   saveas(f5,'Stb3_Kmeans','tif')
                   saveas(f5,'Stb3_Kmeans','pdf')
            elseif ii==3
                   ylabel('Nuclear intensity');                   
                   saveas(f5,'Xbp1_Kmeans')
                   saveas(f5,'Xbp1_Kmeans','tif')
                   saveas(f5,'Xbp1_Kmeans','pdf')
            elseif ii==4
                   ylabel('Nuclear intensity');                   
                   saveas(f5,'Rtg1_Kmeans')
                   saveas(f5,'Rtg1_Kmeans','tif')
                   saveas(f5,'Rtg1_Kmeans','pdf')
            elseif ii==5
                   ylabel('Nuclear intensity');                   
                   saveas(f5,'Gln3_Kmeans')
                   saveas(f5,'Gln3_Kmeans','tif')
                   saveas(f5,'Gln3_Kmeans','pdf')
            elseif ii==6
                   ylabel('Cdc10 index');                   
                   saveas(f5,'Cdc10_Kmeans')
                   saveas(f5,'Cdc10_Kmeans','tif')
                   saveas(f5,'Cdc10_Kmeans','pdf')
            end
            
 %%%%%%%%%%%%----------- obtain heat map ------------------%%%%%%%%%%%%%%%%  
 
groupedMat1=cell(3,1);
for ik=1:Kgroup % loop to concatenate the clusters on top of each other in a single heatmap
pooledMat_noNan = timeseriesM((analysisZ.Splot_infoP(:,9) == ik),:);
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


            if ii==1
                title('Sfp1 nuclear translocation in Cdc10-Clusters')
                   saveas(f6,'Sfp1_heat')
                   saveas(f6,'Sfp1_heat','tif')
                   saveas(f6,'Sfp1_heat','pdf')
            elseif ii==2
                title('Stb3 nuclear translocation in Cdc10-Clusters')
                   saveas(f6,'Stb3_heat')
                   saveas(f6,'Stb3_heat','tif')
                   saveas(f6,'Stb3_heat','pdf')
            elseif ii==3
                title('Xbp1 nuclear translocation in Cdc10-Clusters')
                   saveas(f6,'Xbp1_heat')
                   saveas(f6,'Xbp1_heat','tif')
                   saveas(f6,'Xbp1_heat','pdf')
            elseif ii==4
                title('Rtg1 nuclear translocation in Cdc10-Clusters')
                   saveas(f6,'Rtg1_heat')
                   saveas(f6,'Rtg1_heat','tif')
                   saveas(f6,'Rtg1_heat','pdf')
            elseif ii==5
                title('Gln3 nuclear translocation in Cdc10-Clusters')
                   saveas(f6,'Gln3_heat')
                   saveas(f6,'Gln3_heat','tif')
                    saveas(f6,'Gln3_heat','pdf')
            elseif ii==6
               title('Cdc10-clusters of cells that crossed START at starvation onset')
                   saveas(f6,'Cdc10_heat')
                   saveas(f6,'Cdc10_heat','tif')
                   saveas(f6,'Cdc10_heat','pdf')
            end          
           
end 
   


