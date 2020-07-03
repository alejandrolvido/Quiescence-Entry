

%% phenotypic profiling by machine learning (PPML)

%%%%  This algorithm is a simple way to establish
%%%%  time series clusters, whose main properties (centroids) are stored 
%%%%  and used to find the cluster affiliations of time series from new experiments.
%%%%  See materials and methods Arguello-Miranda, et al. XXXX 

%% define time series period to analyze and normalized load data
clearvars
starti=1;
endti=250;
sr=24;
sr1=5;
%% load data
load('CDC10_timeseries') % normalized times series data
load('Human_labels')% loads manually labelled validation data, which is used to judage how well clustering solutions work
XX=Cdc10_timeseries; % indicates the time series matrix to be clustered
XP = randperm(size(XX,1));% randomizes time series data
Daten=Cdc10_timeseries(XP,:); % randomizes and split time series data
bite=round((size(Daten,1))/10); % randomizes and split time series data
X_learn=Daten((1:(5*bite)),:); % 50% data for learning
% X_val=Daten(((5*bite)+1:8*bite),:); % 30% data for validating, in this case is already labelled in the file "human labels" 
X_test=Daten((8*bite:((10*bite)+1)),:); % 20% data for testing

rng('default') % For reproducibility
X = X_learn;% training set

%% find suggestion for number of clusters using two different methods: silhouette and CalinskiHarabasz.

% by silhouette method
cluster_number1 = evalclusters(X,'kmeans','silhouette','klist',[1:6],'Distance','correlation');

% by CalinskiHarabasz method:
meas=X;
clust = zeros(size(meas,1),6);
for iy=1:6
 clust(:,iy) = kmeans(meas, iy, 'Replicates', 10, 'MaxIter', 10000, 'Display', 'final','Distance','correlation');
end 
cluster_number2 = evalclusters(meas,clust,'CalinskiHarabasz');


display(['optimal cluster number by CalinskiHarabasz method: ' num2str(cluster_number2.OptimalK)]);
display(['optimal cluster number by silhouette method: ' num2str(cluster_number1.OptimalK)]);

%% 

%%%%%%%%% Here biological significance of the number of clusters should be evaluated by biologists before proceding.    


%%

%%this next loop creates the cluster labels by using kmeans on different time windows
%%Correlation is used as distance metric for clustering. 
%%
Kgroup=3; % suggested by silhouette or CalinskiHarabasz methods.  
dista=cell(4,1); % stores types of distance metric used for clustering
dista{1}='correlation';
dista{2}='cityblock';    
dista{3}='sqeuclidean';
dista{4}='cosine';
mmet=numel(dista);
methods=1:mmet;

centroi=cell(mmet,endti); % stores centroid for each clustering interval
mKmns=cell(mmet,endti); % stores the clusters labels for each cells and clustering interval
info_m=cell(2,endti); %  stores the final labels  
MCC=zeros(mmet,endti);
ConfuResult=cell(mmet,(size(X_val,2)-1)); % stores the 

for i=methods

for iw=1:(size(X,2)-1)
[idx,C,sumd,D] = kmeans(X(:,iw:250), Kgroup, 'Replicates', 10, 'MaxIter', 10000, 'Display', 'final','Distance',dista{i});%kmeans
mKmns{iw}=idx;% tabulate(mKmns{i})
centroi{i,iw}=C;
end

%% predicted cluster labels obtained by unsupervised clustering are compared with the manual labels in the validation data set using a multiclass MCC analysis

Idx_predic=zeros(size(X_val));

for ix=1:(size(X_val,2)-2)    
[~,idx] = pdist2(centroi{i,ix}, X_val(:,ix:250),'euclidean','Smallest',1); % prediction of labels using clustering solutions (centroids) for the interval i:end
Idx_predic(:,ix)=idx'; % predicted labels

%% sort mean of each cluster to ensure that cluster labels will be consistent across clustering solutions
info_mn=zeros(size(X_val,1),4);
info_mn(:,1)=(1:size(X_val,1))';
info_mn(:,2)=Idx_predic(:,ix); % 

for ik=1:Kgroup 
pooledMat_noNan = X_val((info_mn(:,2)== ik),:);
% size(pooledMat_noNan )
m1=mean(pooledMat_noNan(:));
ces=find(info_mn(:,2)==ik);
info_mn(ces,3)= m1; 
end

sorces=sort(info_mn(:,3));
unices=unique(sorces)';

for iy=1:3
ces2=find(info_mn(:,3)==unices(1,iy));
info_mn(ces2,4)=iy;
end 

info_m{1,ix}=info_mn(:,4);

%% multiclass MCC to compare predicted vs manual labels
[c_matrix,Result,referenceResult]=confusion.getMatrix(X_info_val,info_m{1,ix});
ConfuResult{i,ix}=Result;
end

%% plot MCC result for the clustering solutions in the interval i:end

for iz=1:(size(X_val,2)-2)
ConfuRe=ConfuResult{i,iz};
MCC(i,iz)=ConfuRe.MatthewsCorrelationCoefficient;
end

if i==max(methods)
plot(MCC')
xlabel('Time point')
ylabel('Multiclass MCC')
legend(dista{1},dista{2},dista{3},dista{4})
title('Matthews Coefficient for clustering solutions')
end
end

%% find best distance metric and best time window for establishing the clusters 
amax=max(MCC(:)); % clustering with the highest Matthew's coefficient is chosen
[a,b]=find(MCC==amax); % chooses the centroids that produced the best solution
C=centroi{a,b};% 

%% visualize and save the optimal centroids to cluster other experiments in the interval i:end that optimizes the clustering

%plot(C') % uncomment to plot
%save('centroids','C'); % uncomment to save

%% used the centroids to cluster new data

X_new = X_test;
[idx_test]=fcentroids(X_new) ;

%% visualize clustering as time series heatmap
timeseriesMat=X_test;% analysisY.SPLOT_6;
X_info_test=idx_test; % asigns cluster affiliations to each time series in analysisY.Splot_info_Cdc10(:,9). 
groupedMat=cell(3,1);
for ik=1:Kgroup % loop to concatenate the clusters on top of each other in a single heatmap
pooledMat_noNan = timeseriesMat((X_info_test== ik),:);
groupedMat{ik}= pooledMat_noNan;
end
 
groupedM_test=[groupedMat{1};groupedMat{2};groupedMat{3}]; % concatenation for heatmap

filteredmap3 = smoothActivityMap(groupedM_test, 'SmoothParam', 0.9, 'UpSample', 1);
figure(4); imagesc(filteredmap3);
colorbar;colormap(jet)
caxis([100 300])
title('Clustered time series')
xlabel('Time (h)');ylabel('Cell Index')
ax = gca;
ax.XTickMode = 'manual';
xticks(starti:sr:endti);
curTick = ax.XTick;
ax.XTickLabel = round((curTick)*sr1/60);


