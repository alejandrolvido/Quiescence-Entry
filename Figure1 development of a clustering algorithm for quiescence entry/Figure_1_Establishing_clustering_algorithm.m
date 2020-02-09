

%%%%  this script is an algorithm to finding the best period to establish
%%%%  clusters for a specific process based on the mattews coefficient between
%%%% clustering assigned by unsupervised algorithms and human-labelled validation data. 

% define period fo the experiment to analyze
clearvars
starti=1;
endti=250;
sr=24;
sr1=5;

load('F1_CDC10_timeseries')
load('Human_labels')% loads manually labelled validation data
XX=Cdc10_timeseries; % indicates the time series matrix to be clustered
XP = randperm(size(XX,1));% randomizes and split time series data
Daten=Cdc10_timeseries(XP,:); % randomizes and split time series data
bite=round((size(Daten,1))/10); % randomizes and split time series data
X_learn=Daten((1:(5*bite)),:); % 50% for learning
% X_val=Daten(((5*bite)+1:8*bite),:); % 30% for validating, previously loaded in this exmaple
X_test=Daten((8*bite:((10*bite)+1)),:); % 20% for testing

rng('default') % For reproducibility
X = X_learn;% training set

Kgroup=3; % can be optimized
centroi=cell(1,250); % keeps centroid for each clustering interval
mKmns=cell(1,250); % keeps the clusters labels for each cells and clustering interval
info_m=cell(2,250); %  keeps the final labels  

%% this loop creates the cluster labels based on kmeans of different state
%% here only one distance metric option (correlation) is displayed for simplicity but notice how changing the distance metric to 
%% sqeuclidean or cityblock will allow to evaluate those options too. 

for i=1:(size(X,2)-1)
[idx,C,sumd,D] = kmeans(X(:,i:250), Kgroup, 'Replicates', 10, 'MaxIter', 10000, 'Display', 'final','Distance','correlation');%kmeans
mKmns{i}=idx;% tabulate(mKmns{i})
centroi{i}=C;
end

% %% optional visualization of centroids
% for i=1:250
% plot(centroi{i}(1,:))
% hold on 
% plot(centroi{i}(2,:))
% hold on 
% plot(centroi{i}(3,:))
% pause
% end

%% predicted cluster labels obtained by unsupervised clustering are compared with the manual labels in the validation data set using a multiclass MCC analysis
ConfuResult=cell(1,(size(X_val,2)-1)); % keeps the 
Idx_predic=zeros(size(X_val));

for i=1:(size(X_val,2)-2)    
[~,idx] = pdist2(centroi{i}, X_val(:,i:250),'euclidean','Smallest',1); % prediction of labels using clustering solutions (centroids) for the interval i:end
Idx_predic(:,i)=idx'; % predicted labels

%% sort mean of each cluster to ensure that cluster labels will be consistent across clustering solutions
info_mn=zeros(size(X_val,1),4);
info_mn(:,1)=(1:size(X_val,1))';
info_mn(:,2)=Idx_predic(:,i); % 

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

info_m{1,i}=info_mn(:,4);

%% multiclass MCC to compare predicted vs manual labels
[c_matrix,Result,referenceResult]=confusion.getMatrix(X_info_val,info_m{1,i});
ConfuResult{i}=Result;
end

%% plot MCC result for the clustering solutions in the interval i:end
MCC=zeros(1,250);
for i=1:(size(X_val,2)-2)
ConfuRe=ConfuResult{i};
MCC(1,i)=ConfuRe.MatthewsCorrelationCoefficient;
end

%% find best initial time point for establishing the clusters and save the corresponding cluster centroids for that interval
plot(MCC)
xlabel('Time point')
ylabel('Multiclass MCC')
title('Matthews Coefficient for clustering intervals')
amax=max(MCC);
besti=find(MCC==amax);
C=centroi{besti};

%% visualize and save the optimal centroids to cluster other experiments in the interval i:end that optimizes the clustering
%plot(C')
%save('centroids','C');

