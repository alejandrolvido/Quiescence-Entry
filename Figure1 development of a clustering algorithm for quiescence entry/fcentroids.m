      function [idx_test]=fcentroids(X_new) % euclidean distance, ignoring NaNs
      load('centroids','')
      [~,idx_test]=pdist2(C,X_new,'euclidean','Smallest',1);