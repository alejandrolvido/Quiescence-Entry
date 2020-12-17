

function []=OAM_180531_plot_subfigure_unvollstandiges4(to_plot_here,colors_here,no_cell_of_interest,current_time,ccell,zoom_in,ccell2,Imagepath,fl_used,fl_used_2,numbM)
              

cell_margin=zoom_in(1,2);
h=figure(1);
set(h,'Position',[5 30 1920 963]);
% set(h,'Position',[2000 -370 1700 1100]);
%set(h,'Position',[2000 -200 1700 1100]);
    
figure(1);subplot(2,2,1);
 
 %to_plot_here=1:260;
for j=1:length(to_plot_here(:,1))
    plot(to_plot_here(j,:),'color',colors_here(j,:));hold on
    plot(current_time,to_plot_here(j,current_time),'Marker','*','LineStyle','none','MarkerSize',20,'LineWidth',4,'Color',colors_here(j,:))
    axis([0 (length(to_plot_here(1,:))) 0 inf]);

end

  
  title(['cell no = ' num2str(no_cell_of_interest)]); 
hold off


%%%%%% Phase Image
% 
ImageName_IP = fullfile( Imagepath ,  sprintf('img_000000%03d_BrightField_000.tif',(current_time)));
IP=imread(ImageName_IP);
% IP=zeros(590,692);
% IP=uint16(IP);
%     imagesc(IP);


%%%%% First Fluorescent Channel

if  fl_used(1,1)==1
ImageName_TFP = fullfile( Imagepath ,  sprintf('img_000000%03d_470 mTFP_000.tif', (current_time)));
I1=imread(ImageName_TFP);
 % I1= ((7*I1)-(medfilt2(I1)));% for stb3
%  I1= 3.5*(2.5*(I1)-(medfilt2(I1))); % for Crz1
 I1=7*(I1-0.7*mean(I1(:)));

else
end

if  fl_used(1,2)==1
ImageName_NG= fullfile( Imagepath ,  sprintf('img_000000%03d_505 mNG_000.tif', (current_time)));%  
%ImageName_NG= fullfile( Imagepath ,  sprintf('img_000000%03d_470 mNG_000.tif', (current_time)));%  

I1=imread(ImageName_NG); 
% I1= ((7*I1)-(medfilt2(I1)));
I1=7*(I1-0.7*mean(I1(:)));
else
end

if  fl_used(1,3)==1
ImageName_K= fullfile( Imagepath ,  sprintf('img_000000%03d_555 mKok_000.tif', (current_time)));%  
I1=imread(ImageName_K);
%  I1=7*(I1-0.7*mean(I1(:)));%sfp1
 I1=10*(I1-0.4*mean(I1(:)));%gln3
 %I1= ((7*I1)-(medfilt2(I1)));
%I1= 6.7*((1.7*I1)-(medfilt2(I1)));
else
end

if  fl_used(1,4)==1
ImageName_Ruby= fullfile( Imagepath ,  sprintf('img_000000%03d_555 mRuby3_000.tif', (current_time)));%  
I1=imread(ImageName_Ruby);
% I1= 4*(10*(I1)-(medfilt2(I1)));%Msn2
 I1=7*(I1-0.7*mean(I1(:)));%mig1
%  I1=I1-0.95*I1;
else
end

if  fl_used(1,5)==1
ImageName_N= fullfile( Imagepath ,  sprintf('img_000000%03d_615 nm_000.tif', (current_time)));%  
I1=imread(ImageName_N);
%I1= ((3*I1)-(medfilt2(I1)));% for stb3
I1= ((1.8*I1));%-(4*medfilt2(I1)));
else
end


if  fl_used(1,6)==1
ImageName_C= fullfile( Imagepath ,  sprintf('img_000000%03d_470 GFP_000.tif', (current_time)));%  
I1=imread(ImageName_C);
I1=7*(I1-0.7*mean(I1(:)));
%  I1=zeros(590,692);
%  I1=uint16(I1);
else
end

if  fl_used(1,7)==1
ImageName_C= fullfile( Imagepath ,  sprintf('img_000000%03d_470 newGFP_000.tif', (current_time)));%  
I1=imread(ImageName_C);
I1=7*(I1-0.7*mean(I1(:)));
else
end

%%%%%%% Second fluorescent channel

if  fl_used_2(1,1)==1
ImageName_TFP = fullfile( Imagepath ,  sprintf('img_000000%03d_470 mTFP_000.tif', (current_time)));
I2=imread(ImageName_TFP);
I2= ((7*I2)-(medfilt2(I2)));% for stb3
% %I2= 1.5*(4*(I2)-(medfilt2(I2))); % for Crz1
% I2=zeros(590,692);
% I2=uint16(I2);
else
end

if  fl_used_2(1,2)==1
ImageName_NG= fullfile( Imagepath ,  sprintf('img_000000%03d_505 mNG_000.tif', (current_time)));%  
I2=imread(ImageName_NG); 
%I2=0.5*I2;
I2= 0.5*(2*(I2)-(medfilt2(I2))); % for Crz1

else
end

if  fl_used_2(1,3)==1
ImageName_K= fullfile( Imagepath ,  sprintf('img_000000%03d_555 mKok_000.tif', (current_time)));%  
I2=imread(ImageName_K);
I2= ((7*I2)-(medfilt2(I2)));
else
end

if  fl_used_2(1,4)==1
ImageName_Ruby= fullfile( Imagepath ,  sprintf('img_000000%03d_555 mRuby3_000.tif', (current_time)));%  
I2=imread(ImageName_Ruby);
% I2= 15*(10*(I2)-(medfilt2(I2)));%Msn2
 I2=10*(I2-0.9*mean(I2(:)));%mig1
else
end

if  fl_used_2(1,5)==1
ImageName_N= fullfile( Imagepath ,  sprintf('img_000000%03d_615 nm_000.tif', (current_time)));%  
I2=imread(ImageName_N);
else
end

if  fl_used_2(1,6)==1
ImageName_C= fullfile( Imagepath ,  sprintf('img_000000%03d_470 GFP_000.tif', (current_time)));%  
I2=imread(ImageName_C);
I2=2*(I2-0.7*mean(I2(:)));

else
end


if  fl_used_2(1,7)==1
ImageName_C= fullfile( Imagepath ,  sprintf('img_000000%03d_470 newGFP_000.tif', (current_time)));%  
I2=imread(ImageName_C);
I2=7*(I2-0.7*mean(I2(:)));
else
end

%%%%%%%%%% subplotting
if zoom_in(1,1)==1
    
  [x_cn,y_cn]=get_wind_coord1(ccell,cell_margin);
   outline=48.*uint16(bwmorph(ccell2,'remove',inf))+max(I1(:)).*uint16(bwmorph(ccell,'remove'));
   outline2=48.*uint16(bwmorph(ccell2,'remove',inf))+max(IP(:)).*uint16(bwmorph(ccell,'remove'));
   outline3=48.*uint16(bwmorph(ccell2,'remove',inf))+max(I2(:)).*uint16(bwmorph(ccell,'remove'));
   outline_I1=(outline(y_cn,x_cn));
   outlineIP=(outline2(y_cn,x_cn));
   outline_I2=(outline3(y_cn,x_cn));
       
   subplot(2,2,2);
   imagesc(0.5*outline_I1+10*I1(y_cn,x_cn));title('FL 1');
   text(5,5,['time = ' num2str((current_time))],'color',[1 1 0.5]);
   
   subplot(2,2,3);
   imagesc(outlineIP+IP(y_cn,x_cn));title('phase');
   text(10,10,['time = ' num2str((current_time))],'color',[0 1 0]);
     
   subplot(2,2,4);
   imagesc(0.5*outline_I2+10*I2(y_cn,x_cn));title('FL 2');
   text(10,10,['time = ' num2str((current_time))],'color',[1 1 0]);
    
else

    
end
end