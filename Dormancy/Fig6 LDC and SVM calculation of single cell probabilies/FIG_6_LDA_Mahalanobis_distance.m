%%% this code uses an optimized LDA for each timepoint and measures the
%%% cluster separation by measuring the mahalanobis distance to one of the
%%% clusters

clearvars

%% define channels
load('F5_SVM_Cdc10_time series')

nn=8; % column where cluster labels are stored
nn1=7;% column where the number of biological replicate is stored
predic1=Cdc10_time_series; % matrix corresponding to the channel to to predict

maja1=zeros(size(predic1)); %  column that stores mahalanobis of cluster 1 distance to cluster 2
maja2=zeros(size(predic1)); % column that stores mahalanobis of cluster 2 distance to cluster 2
Mdl1=cell(1,259); % stores the LDA model for ech time point

%% loop to create a linear discriminant classifier for each time point and calculate Mahalanobis distance to cluster 2

%%%% notice this step might take long when all time points are analized 
for i=1:259
rng(1)
Mdl = fitcdiscr(predic1(:,i),Cdc10_time_series_info(:,8),'OptimizeHyperparameters','auto',...
    'HyperparameterOptimizationOptions',...
    struct('Holdout',0.3,'AcquisitionFunctionName','expected-improvement-plus'));
maha= mahal(Mdl,predic1(:,i));
maja1(:,i)=maha(:,1);
maja2(:,i)=maha(:,2);
Mdl1{i}=Mdl;
close all
end

Acells=find(Cdc10_time_series_info(:,nn)==2); % find cells in cluster 1
Bcells=find(Cdc10_time_series_info(:,nn)==3); % find cells in cluster 2

%% visualise as heatmap 
f1= figure; % 
figtmp = imagesc([maja1(Bcells,:); maja1(Acells,:)]);%
colorbar;colormap(jet)
caxis([0 20])

    titlex=('Mahalanobis distance to cluster 2');         
    title(titlex)
    xlabel('Time point');
    ylabel('Cell index');

%% visualize as average with 95% confidence intervals
numCells = 1:length(predic1);
chamId = Cdc10_time_series_info(:,nn1); 
CH_no=max(chamId);
meanArr_perCham = cell(2,1);
timeseriesM=maja1;
Kgroup=1:2;

for i=1:2 % loop to obtain the average of each cluster per lane 
    meanArr_perCham{i} = nan(length(Kgroup), length(timeseriesM(1,:)));
    for j=1:CH_no
        ind = (((Cdc10_time_series_info(:,nn)-1) == i) & (chamId == j));
        tmp = timeseriesM(ind,:);
        tmp1 = nanmean(tmp, 1);
        meanArr_perCham{i}(j, :) = tmp1;
    end
end

meanOfChambers = nan(length(Kgroup), length(meanArr_perCham{1}));
for i=Kgroup % loop to obtain the mean of each lane, biological replicate
    meanOfChambers(i,:) = nanmean(meanArr_perCham{i}, 1);
end
t=1:length(meanArr_perCham{1});

colrs=['b'; 'r'; 'k'; 'g'; 'c';'m';'y'];
limX=[24 259];
sr1=24;
sr2=5;
starti=1;
endti=259;

f5=figure;
group=1:2;
for ik4=group  % loop to plot the average curve for each cluster
    sem = std(meanArr_perCham{ik4}, [], 1, 'omitnan')./sqrt(size(meanArr_perCham{ik4},1));
    hold on
    s1 = shadedErrorBarV2(t, meanOfChambers(ik4,:), 2*sem, 'lineprops', colrs(ik4));
    hold on
end
              
    titlex=('Separation of clusters by Cdc10 (Mahalanobis distance)');         
    title(titlex)
    xlim(limX)
    xlabel('Time point');
    ylabel('Mahalanobis Distance to cluster 2');  

%% end




