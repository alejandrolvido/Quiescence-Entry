function [Mat_compress,GR_Matrix] = find_maxgrowthrate(Mat1, nfactor)

Mat_compress=zeros(size(Mat1));

%%%%outliers detection and imputation
for i=1:length(Mat1(:,1))
MN=nanmean(Mat1(i,:));
ST=nanstd(Mat1(i,:));
OUT=find(Mat1(i,:)>MN+3*ST);
if OUT~=0
for j=OUT
   if j==1
Mat1(i,j)=Mat1(i,j+1);
   elseif j==max(length(Mat1(1,:)))
Mat1(i,j)=Mat1(i,j-1);
   else
Mat1(i,j)=(Mat1(i,j-1)+Mat1(i,j+1))/2;
   end
end
else
end
end
%%%%% replace NAN by 0 for compression
A=isnan(Mat1);
Mat1(A)=0;

if nfactor==2  
    G2_matrix1=zeros(size(Mat1,1), (size(Mat1,2)-1)/2);
for j=1:length(Mat1(:,1))
    b=1;
    for i=1:2:length(Mat1(1,:))-1             

        if (Mat1(j,i))==0
            G2_matrix1(j,b)=(Mat1(j,i+1));
        elseif Mat1(j,i+1)==0
            G2_matrix1(j,b)=(Mat1(j,i));
        else
            G2_matrix1(j,b)=((Mat1(j,i))+(Mat1(j,i+1)))/2;
        end
        b=b+1;
    end
end

%%%%% replace NAN by 0 for compression
A=isnan(G2_matrix1);
G2_matrix1(A)=0;
Mat_compress=G2_matrix1;

elseif nfactor==4
    %%%%% compression I
G2_matrix1=zeros(size(Mat1,1), (size(Mat1,2)-1)/2);
for j=1:length(Mat1(:,1))
    b=1;
    for i=1:2:length(Mat1(1,:))-1             

        if (Mat1(j,i))==0
            G2_matrix1(j,b)=(Mat1(j,i+1));
        elseif Mat1(j,i+1)==0
            G2_matrix1(j,b)=(Mat1(j,i));
        else
            G2_matrix1(j,b)=((Mat1(j,i))+(Mat1(j,i+1)))/2;
        end
        b=b+1;
    end
end

%%%%% replace NAN by 0 for compression
A=isnan(G2_matrix1);
G2_matrix1(A)=0;

%%%%% compression II
G3_matrix1=zeros(size(G2_matrix1,1), (size(G2_matrix1,2)-1)/2);
for j=1:length(G2_matrix1(:,1))
    b=1;
    for i=1:2:length(G2_matrix1(1,:))-1             

        if (G2_matrix1(j,i))==0
            G3_matrix1(j,b)=(G2_matrix1(j,i+1));
        elseif G2_matrix1(j,i+1)==0
            G3_matrix1(j,b)=(G2_matrix1(j,i));
        else
            G3_matrix1(j,b)=((G2_matrix1(j,i))+(G2_matrix1(j,i+1)))/2;
        end
        b=b+1;
    end
end

A=isnan(G3_matrix1);
G3_matrix1(A)=0;

Mat_compress=G2_matrix1;

elseif nfactor==8
    
      %%%%% compression I
G2_matrix1=zeros(size(Mat1,1), (size(Mat1,2)-1)/2);
for j=1:length(Mat1(:,1))
    b=1;
    for i=1:2:length(Mat1(1,:))-1             

        if (Mat1(j,i))==0
            G2_matrix1(j,b)=(Mat1(j,i+1));
        elseif Mat1(j,i+1)==0
            G2_matrix1(j,b)=(Mat1(j,i));
        else
            G2_matrix1(j,b)=((Mat1(j,i))+(Mat1(j,i+1)))/2;
        end
        b=b+1;
    end
end

A=isnan(G2_matrix1);
G2_matrix1(A)=0;

%%%%% compression II
G3_matrix1=zeros(size(G2_matrix1,1), (size(G2_matrix1,2)-1)/2);
for j=1:length(G2_matrix1(:,1))
    b=1;
    for i=1:2:length(G2_matrix1(1,:))-1             

        if (G2_matrix1(j,i))==0
            G3_matrix1(j,b)=(G2_matrix1(j,i+1));
        elseif G2_matrix1(j,i+1)==0
            G3_matrix1(j,b)=(G2_matrix1(j,i));
        else
            G3_matrix1(j,b)=((G2_matrix1(j,i))+(G2_matrix1(j,i+1)))/2;
        end
        b=b+1;
    end
end

A=isnan(G3_matrix1);
G3_matrix1(A)=0;

%%%%% compression III
G4_matrix1=zeros(size(G3_matrix1,1), (size(G3_matrix1,2)-1)/2);
for j=1:length(G3_matrix1(:,1))
    b=1;
    for i=1:2:length(G3_matrix1(1,:))-1             

        if (G3_matrix1(j,i))==0
            G4_matrix1(j,b)=(G3_matrix1(j,i+1));
        elseif G3_matrix1(j,i+1)==0
            G4_matrix1(j,b)=(G3_matrix1(j,i));
        else
            G4_matrix1(j,b)=((G3_matrix1(j,i))+(G3_matrix1(j,i+1)))/2;
        end
        b=b+1;
    end
end

A=isnan(G4_matrix1);
G4_matrix1(A)=0;

Mat_compress=G4_matrix1;

elseif nfactor==16
          %%%%% compression I
G2_matrix1=zeros(size(Mat1,1), (size(Mat1,2)-1)/2);
for j=1:length(Mat1(:,1))
    b=1;
    for i=1:2:length(Mat1(1,:))-1             

        if (Mat1(j,i))==0
            G2_matrix1(j,b)=(Mat1(j,i+1));
        elseif Mat1(j,i+1)==0
            G2_matrix1(j,b)=(Mat1(j,i));
        else
            G2_matrix1(j,b)=((Mat1(j,i))+(Mat1(j,i+1)))/2;
        end
        b=b+1;
    end
end

A=isnan(G2_matrix1);
G2_matrix1(A)=0;

%%%%% compression II
G3_matrix1=zeros(size(G2_matrix1,1), (size(G2_matrix1,2)-1)/2);
for j=1:length(G2_matrix1(:,1))
    b=1;
    for i=1:2:length(G2_matrix1(1,:))-1             

        if (G2_matrix1(j,i))==0
            G3_matrix1(j,b)=(G2_matrix1(j,i+1));
        elseif G2_matrix1(j,i+1)==0
            G3_matrix1(j,b)=(G2_matrix1(j,i));
        else
            G3_matrix1(j,b)=((G2_matrix1(j,i))+(G2_matrix1(j,i+1)))/2;
        end
        b=b+1;
    end
end

A=isnan(G3_matrix1);
G3_matrix1(A)=0;

%%%%% compression III
G4_matrix1=zeros(size(G3_matrix1,1), (size(G3_matrix1,2)-1)/2);
for j=1:length(G3_matrix1(:,1))
    b=1;
    for i=1:2:length(G3_matrix1(1,:))-1             

        if (G3_matrix1(j,i))==0
            G4_matrix1(j,b)=(G3_matrix1(j,i+1));
        elseif G3_matrix1(j,i+1)==0
            G4_matrix1(j,b)=(G3_matrix1(j,i));
        else
            G4_matrix1(j,b)=((G3_matrix1(j,i))+(G3_matrix1(j,i+1)))/2;
        end
        b=b+1;
    end
end

A=isnan(G4_matrix1);
G4_matrix1(A)=0;
    
%%%%% compression IV
G5_matrix1=zeros(size(G4_matrix1,1), (size(G4_matrix1,2))/2);
for j=1:length(G4_matrix1(:,1))
    b=1;
    for i=1:2:length(G4_matrix1(1,:))-1             

        if (G4_matrix1(j,i))==0
            G5_matrix1(j,b)=(G4_matrix1(j,i+1));
        elseif G4_matrix1(j,i+1)==0
            G5_matrix1(j,b)=(G4_matrix1(j,i));
        else
            G5_matrix1(j,b)=((G4_matrix1(j,i))+(G4_matrix1(j,i+1)))/2;
        end
        b=b+1;
    end
end

A=isnan(G5_matrix1);
G5_matrix1(A)=0;
    
Mat_compress=G5_matrix1;
end

V_mat=Mat_compress;
GR_Matrix=zeros(size(Mat_compress));
for i=1:length(V_mat(:,1))
       for j=2:length(V_mat(1,:))
    GR_Matrix(i,j-1)=(V_mat(i,j)-V_mat(i,j-1));
       end
end       
       infCor=GR_Matrix==inf;
        GR_Matrix(infCor)=NaN;
%         
%         GR_Matrix=zeros(size(GR_Matrix));
% %%% Outlixers II ixmputatixons
% for ix=1:length(GR_Matrix1(:,1))
% MN=nanmean(GR_Matrix1(ix,:));
% ST=nanstd(GR_Matrix1(ix,:));
% OUT=find(GR_Matrix1(ix,:)>MN+3*ST | GR_Matrix1(ix,:)<MN-3*ST);
% if OUT~=0
% for j=OUT
%    if j==1
% GR_Matrix(ix,j)=GR_Matrix1(ix,j+1);
%    elseif j==max(length(GR_Matrix1(1,:)))
% GR_Matrix(ix,j)=GR_Matrix1(ix,j-1);
%    else
% GR_Matrix(ix,j)=(GR_Matrix1(ix,j-1)+GR_Matrix1(ix,j+1))/2;
%    end
% end
% else
% end        
end