
%%%% this code exemplies how to calculate single cell probabilities based
%%%% using a SVM trained on the final state of a time series
%%%% dotted line show pior probabily treshold 

clearvars

%% define channels
load('F5_SVM_Cdc10_time series')

predic=cell(1,259);
tp=259;
colrs=['b'; 'r'; 'k'; 'g'; 'c';'m';'y'];
limX=[24 259];
limY=[0 1];
sr1=24;
sr2=5;
starti=24; %
endti=259; %

% Channel used for calcualting single cell probabilities

 predic1=Cdc10_time_series(:,259); % final timepoints used to train the SVM 
 SVMModel = fitcsvm(predic1,Cdc10_time_series_info(:,8),'KernelFunction','rbf',...
    'Standardize',true,'ClassNames',{'2','3'});

ccellL=zeros(size(Cdc10_time_series));
ccellL2=zeros(size(Cdc10_time_series));
%% loop to obtain posterior probabilities for each time point based on the trained SVM
 ccell=cell(1,259);
for i=1:259
ScoreSVMModel = fitPosterior(SVMModel,Cdc10_time_series(:,i),Cdc10_time_series_info(:,8));
[label,score] = predict(ScoreSVMModel,Cdc10_time_series(:,i));
ccell{i}=score;
end

%% loop to gather posterior probabilities
for i=1:259
     ccellL(:,i)=ccell{i}(:,1);
    ccellL2(:,i)=ccell{i}(:,2);
end

%% plotting average posterior probabilities per biological replicate with 95% confidence intervals 
numCells = 1:length(Cdc10_time_series(:,i));
chamId = Cdc10_time_series_info(:,7); 
CH_no=max(chamId);
 meanArr_perCham = cell(2,1);
timeseriesM=ccellL;
Kgroup=1:2;

for i=1:2 % loop to obtain the average of each cluster per lane in the microfluidic chamber (replicate)
    meanArr_perCham{i} = nan(length(Kgroup), length(timeseriesM(1,:)));
    for j=1:CH_no
        ind = (((Cdc10_time_series_info(:,8)-1) == i) & (chamId == j));
        tmp = timeseriesM(ind,:);
        tmp1 = nanmean(tmp, 1);
        meanArr_perCham{i}(j, :) = tmp1;
    end
end

meanOfChambers = nan(length(Kgroup), length(meanArr_perCham{1}));
for i=Kgroup % loop to obtain the mean of each lane
    meanOfChambers(i,:) = nanmean(meanArr_perCham{i}, 1);
end
t=1:length(meanArr_perCham{1});

f5=figure;
group=1:2;
for ik4=group  % loop to plot the average curve for each cluster
    sem = std(meanArr_perCham{ik4}, [], 1, 'omitnan')./sqrt(size(meanArr_perCham{ik4},1));
    hold on
    s1 = shadedErrorBarV2(t, meanOfChambers(ik4,:), 2*sem, 'lineprops', colrs(ik4));
    hold on
end
              
    titlex=('Average Probability of division');         
    title(titlex)
    xlim(limX)
    ylim(limY)
    xlabel('Time (h)');
    ax = gca;
    ax.XTickMode = 'manual';
    xticks((starti:sr1:endti));
    curTick = ax.XTick;
    ax.XTickLabel = round((curTick)*sr2/60)-starti*sr2/60;
    line(1:259,repelem(0.8737,259),'Color','red','LineStyle','--');
saveas(f5,titlex)
saveas(f5,titlex,'pdf')

timeseriesM=ccellL2;
for i=1:2 % loop to obtain the average of each cluster per lane 
    meanArr_perCham{i} = nan(length(Kgroup), length(timeseriesM(1,:)));
    for j=1:CH_no
        ind = (((Cdc10_time_series_info(:,8)-1) == i) & (chamId == j));
        tmp = timeseriesM(ind,:);
        tmp1 = nanmean(tmp, 1);
        meanArr_perCham{i}(j, :) = tmp1;
    end
end

meanOfChambers = nan(length(Kgroup), length(meanArr_perCham{1}));
for i=Kgroup % loop to obtain the mean of each lane
    meanOfChambers(i,:) = nanmean(meanArr_perCham{i}, 1);
end
t=1:length(meanArr_perCham{1});

f6=figure;
group=1:2;
for ik4=group  % loop to plot the average curve for each cluster
    sem = std(meanArr_perCham{ik4}, [], 1, 'omitnan')./sqrt(size(meanArr_perCham{ik4},1));
    hold on
    s1 = shadedErrorBarV2(t, meanOfChambers(ik4,:), 2*sem, 'lineprops', colrs(ik4));
    hold on
end
             
    titlex=('average Probability of cell cycle arrest');                 
    title(titlex)
    xlim(limX)
    ylim(limY)
    xlabel('Time (h)');
    ax = gca;
    ax.XTickMode = 'manual';
    xticks(starti:sr1:endti);
    curTick = ax.XTick;
    ax.XTickLabel = round((curTick)*sr2/60)-starti*sr2/60;
    line(1:259,repelem((1-0.8737),259),'Color','red','LineStyle','--');
saveas(f6,titlex)
saveas(f6,titlex,'pdf')


