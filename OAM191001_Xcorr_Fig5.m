
%% Cross correlation analysis code for single cells time series


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
load('Toy_data_Fig5')

MatVecX=analysis.SPLOT_1; % primary matrix to use as reference for cross-correlations with all other channels, i.e. Sfp1. 
colrs=['r'; 'b'; 'k'; 'g'; 'c';'m';'y']; % colors for plotting
Kgroup=5; % number of clusters
markers=(2); % secondary matrix to compute the cross-correlations, 1=sfp1, 2=stb3,3=Xbp1,4=Rtg1,5=Gln3,6=Cdc10 
savef=0; % saving if 1
t=-33:1:33; % sliding windown time used for cross-correlations, usually 25%-33% of the whole experiment duration. 
sr=5; % sampling rate
Xzero=34; % set zero in plots
limY=[-0.8 0.8]; % defines y axis according to marker
caxi=([-0.4 0.4]); % defines color bar in heatmap 
path_save='your path to save'; % save 
for i=markers % loop through all the channels that will be crosscorrelated to the abovedefined reference matrix 

% declare the matrices that will be crosscorrelated to the reference
% matrix 

 if i==1
    MatVecY=analysis.SPLOT_1;
    titel1='Sfp1 Xcorr heat map';
    titel2='Population average Sfp1 Xcorr';
    titel3='Average Sfp1 Xcorr by k-cluster';
 end

 if i==2
    MatVecY=analysis.SPLOT_2;
    titel1='Stb3 Xcorr heat map';
    titel2='Population average Stb3  Xcorr';
    titel3='average Stb3  Xcorr by k-cluster';
 end

 if i==3
    MatVecY=analysis.SPLOT_3;
    titel1='Xbp1 Xcorr heat map';
    titel2='Population average Xbp1  Xcorr';
    titel3='average Xbp1  Xcorr by k-cluster';
 end

 if i==4
    MatVecY=analysis.SPLOT_4;
    titel1='Rtg1 Xcorr heat map';
    titel2='Population average Rtg1  Xcorr';
    titel3='average Rtg1  Xcorr by k-cluster';
 end

 if i==5
    MatVecY=analysis.SPLOT_5;
    titel1='Gln3 Xcorr heat map';
    titel2='Population average Gln3 Xcorr';
    titel3='average Gln3 Xcorr by k-cluster';
 end

 if i==6
    MatVecY=analysis.SPLOT_6;
    titel1='Cdc10 Xcorr heat map';
    titel2='Population average Cdc10  Xcorr';
    titel3='average Cdc10  Xcorr by k-cluster';
 end
 

% calculates Xcorr 

xCorrMap = nanXcorrMaps(MatVecX, MatVecY);

% loop calculating the Xcorr in each K-cluster

Ch_NO=max(analysis.Splot_info(:,7));
groupedMat=cell(Kgroup,1);
for ik=1:Kgroup
pooled = xCorrMap(analysis.Splot_info(:,8)==ik,:);
groupedMat{ik} = pooled;
end
groupedM=[groupedMat{1};groupedMat{2};groupedMat{3};groupedMat{4};groupedMat{5}]; % simple concatenation to allow for custom stacking of clusters 

% visualize as heatmap
inputmap1 = groupedM;
filteredmap = smoothActivityMap(inputmap1, 'SmoothParam', 0.9, 'UpSample', 1);
f4 = figure;
figtmp = imagesc(filteredmap);
colorbar;colormap(jet);
caxis(caxi)
figtmp.AlphaData = 1-isnan(inputmap1);
title(titel1)
xticks( [Xzero-24  Xzero-12 Xzero Xzero+12 Xzero+24]);
xlabel('Time (h)');ylabel('Xcorr')
ax = gca; 
xlim([1 67])   
ax.XTickMode = 'manual';
ax.XTickLabel=( [-2 -1 0 1 2]);

% save if 1
 if savef==1
export_path=[path_save titel1];
hgexport(gcf, export_path, hgexport('factorystyle'), 'Format', 'pdf');
 else
 end

 
%% loop calculating the Xcorr per cluster in each lane (replicate)
numCells = size(xCorrMap, 1);
Lane_Id = analysis.Splot_info(:,7); 

meanArr_perCham = cell(Kgroup,1);
for iyy=1:Ch_NO
    meanArr_perCham{iyy} = nan(5, size(xCorrMap, 2));
    
    for j=1:Kgroup
        ind = ((analysis.Splot_info(:,8) == j) & (Lane_Id == iyy));
        tmp = xCorrMap(ind,:);
        tmp1 = nanmean(tmp, 1);
        meanArr_perCham{iyy}(j, :) = tmp1;
    end
end

%% plot population mean curve 

% calculates mean of each lane
meanOfLanes = nan(Kgroup, length(meanArr_perCham{1}));
for ik2=1:Kgroup
    meanOfLanes(ik2,:) = nanmean(meanArr_perCham{ik2}, 1);
end

%%
f5=figure ;

    sem = std(meanOfLanes, [], 1, 'omitnan')./sqrt(size(meanArr_perCham,1));
    hold on
    s1 = shadedErrorBarV2(t, nanmean(meanOfLanes(:,:)), 2*sem, 'lineprops', colrs(1));
  
title(titel2)
xlabel('Time (h)');ylabel('Xcorr')
ax = gca;
 xlim([-33 33])   
ax.XTickMode = 'manual';
 xticks( [-24 -12 0 12 24]);
curTick = ax.XTick;
ax.XTickLabel = ((curTick)*sr/60);%round
xlabel('Time (h)');
ylim(limY);
 if savef==1
%export_path=['A:\scores\Orlando\Cyclic\Final Analysis with GR\' titel2];
export_path=[path_save titel2];
hgexport(gcf, export_path, hgexport('factorystyle'), 'Format', 'pdf');
 else
 end 
    
%% plot mean curves per cluster
f6=figure ;

for ik4=1:Kgroup
    sem = std(meanArr_perCham{ik4}, [], 1, 'omitnan')./sqrt(size(meanArr_perCham{ik4},1));
    hold on
    s1 = shadedErrorBarV2(t, meanOfLanes(ik4,:), 2*sem, 'lineprops', colrs(ik4));
    hold on
end
              
title(titel3)
xlabel('Time (h)');ylabel('Xcorr')
ax = gca;
 xlim([-33 33])   
ax.XTickMode = 'manual';
 xticks( [-24 -12 0 12 24]);
curTick = ax.XTick;
ax.XTickLabel = ((curTick)*sr/60);%round
xlabel('Time (h)');
ylim(limY);

 if savef==1
export_path=[path_save titel3];
hgexport(gcf, export_path, hgexport('factorystyle'), 'Format', 'pdf');
 else
 end     

end


