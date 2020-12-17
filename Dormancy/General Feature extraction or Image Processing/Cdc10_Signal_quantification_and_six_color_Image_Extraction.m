%%% this script is used to extract the Cdc10 Signal and 
%%% the nuclear intensity of stress-transcription factors in up to six fluorescent channels plus 
%%% phase contrast 

clearvars
parpool('local',12)

type='.tif';
path_h='your path';
exp_name='your experiment';
positions=[0];% Field of view
disp(positions)

fl_used=[1 1 1 1 1 1];% channels for analysis, in order of appearance TFP GFP/mNG mKok mRuby mNeptune CYOFP
% =========================================================================
bin_used=1;
IT_mean_modifier= 1;
G_mean_modifier= 1;
K_mean_modifier= 1;
R_mean_modifier= 1;
FR_mean_modifier= 1;
peak_cutoff=0.75;
final_tp=30;



for pos_here=positions
    ptmp_h=pos_here;
    path_seg=path_h;
    load([path_seg exp_name '/' exp_name '_pos' num2str(ptmp_h)])
    
    Imagepath=([path_h exp_name '/Pos' num2str(pos) '/']); % path to images          
    
    % define image size in pixels 
    x_size =520;   %1040;
    y_size =692;   %1388;
    numbM=final_tp;
    appr_vol                      =zeros(no_obj,numbM);
    
    %---------------------------------------------------------------------
  % Allocation of features extracted from all fluorophores regardless of the
  % protein tagged
    
    %%%%%%%--------------- mTFP1 feature allocation ----------------%%%%%%%
    
    if fl_used(1,1)==1
        all_obj.med_T                 =zeros(1,numbM);
        mean_IT_int_per_area_T        =zeros(no_obj,numbM);
        IT_Conc_T                     =zeros(no_obj,numbM);
        max_nucl_int_IT               =zeros(no_obj,numbM);
    
        %%%-------------------- features extracted for mTFP --------------------%%%

        
    nucl_area_TFP                  =zeros(no_obj,numbM);
    cyt_area_TFP                   =zeros(no_obj,numbM);
    mean_TFP_int_per_area_C        =zeros(no_obj,numbM);
    mean_TFP_int_per_area_N        =zeros(no_obj,numbM);
    TFP_Conc_T                     =zeros(no_obj,numbM);
    TFP_Conc_C                     =zeros(no_obj,numbM);
    TFP_Conc_N                     =zeros(no_obj,numbM);
    TFP_mean_int_N                 =zeros(no_obj,numbM);
    TFP_mean_int_C                 =zeros(no_obj,numbM);
    tot_TFP_fl_cyt                 = zeros(no_obj,numbM);
    tot_TFP_fl_nuc                 = zeros(no_obj,numbM);
        
    end
    
    %%%%%%%------------- mNeonGreen feature allocation -------------%%%%%%%
    
    if fl_used(1,2)==1
        all_obj.med_G                 =zeros(1,numbM);
        mean_Green_int_per_area_T     =zeros(no_obj,numbM);
        Green_Conc_T                  =zeros(no_obj,numbM);
        max_nucl_int_Green            =zeros(no_obj,numbM);
        
        
%%%----------------- features extracted for NeonGreen ------------------%%%
nucl_area_Neon_Green                  =zeros(no_obj,numbM);    
cyt_area_Neon_Green                  =zeros(no_obj,numbM);
    mean_Neon_Green_int_per_area_C        =zeros(no_obj,numbM);
    mean_Neon_Green_int_per_area_N        =zeros(no_obj,numbM);
    Neon_Green_Conc_T                     =zeros(no_obj,numbM);
    Neon_Green_Conc_C                     =zeros(no_obj,numbM);
    Neon_Green_Conc_N                     =zeros(no_obj,numbM);
    nucl_Neon_Green                      =false(x_size,y_size,numbM);
    Neon_Green_mean_int_N               =zeros(no_obj,numbM);
    Neon_Green_mean_int_C               =zeros(no_obj,numbM);
    tot_Neon_Green_fl_cyt               = zeros(no_obj,numbM);
    tot_Neon_Green_fl_nuc               = zeros(no_obj,numbM);
        
    end
    
    %%%%%%%---------------- mKOk feature allocation ----------------%%%%%%%
    
    if fl_used(1,3)==1
        all_obj.med_K                 =zeros(1,numbM);
        mean_K_int_per_area_T         =zeros(no_obj,numbM);
        K_Conc_T                      =zeros(no_obj,numbM);
        max_nucl_int_K                =zeros(no_obj,numbM);
  
 %%%-------------------- features extracted for mKOk --------------------%%%
    nucl_area_KOK                  =zeros(no_obj,numbM);
    cyt_area_KOK                   =zeros(no_obj,numbM);
    mean_KOK_int_per_area_C        =zeros(no_obj,numbM);
    mean_KOK_int_per_area_N        =zeros(no_obj,numbM);
    KOK_Conc_T                     =zeros(no_obj,numbM);
    KOK_Conc_C                     =zeros(no_obj,numbM);
    KOK_Conc_N                     =zeros(no_obj,numbM);
    KOK_mean_int_N                 =zeros(no_obj,numbM);
    KOK_mean_int_C                 =zeros(no_obj,numbM);
    tot_KOK_fl_cyt                 = zeros(no_obj,numbM);
    tot_KOK_fl_nuc                 = zeros(no_obj,numbM);
        
        
    end
    
    %%%%%%%------ mRuby3 or mScarlet-I feature allocation ----------%%%%%%%
    
    if fl_used(1,4)==1
        all_obj.med_mRu               =zeros(1,numbM);
        mean_Ru_int_per_area_T        =zeros(no_obj,numbM);
        Ru_Conc_T                     =zeros(no_obj,numbM);
        max_nucl_int_Ru               =zeros(no_obj,numbM);   
        
        %%%------------ features extracted for mRuby or mScarlet-I ------------%%%
    nucl_area_Ruby                  =zeros(no_obj,numbM);
    cyt_area_Ruby                   =zeros(no_obj,numbM);
    mean_Ruby_int_per_area_C        =zeros(no_obj,numbM);
    mean_Ruby_int_per_area_N        =zeros(no_obj,numbM);
    Ruby_Conc_T                     =zeros(no_obj,numbM);
    Ruby_Conc_C                     =zeros(no_obj,numbM);
    Ruby_Conc_N                     =zeros(no_obj,numbM);
    Ruby_mean_int_N                 =zeros(no_obj,numbM);
    Ruby_mean_int_C                 =zeros(no_obj,numbM);
    tot_Ruby_fl_cyt                 = zeros(no_obj,numbM);
    tot_Ruby_fl_nuc                 = zeros(no_obj,numbM);
    end
    
    %%%%%%%------------ mNeptune2.5 feature allocation -------------%%%%%%%
    
    if fl_used(1,5)==1
        all_obj.med_mCa               =zeros(1,numbM);
        mean_Nep_int_per_area_T       =zeros(no_obj,numbM);
        Nept_Conc_T                   =zeros(no_obj,numbM);
        max_nucl_int_Nep              =zeros(no_obj,numbM);
        
        
        %%%------------------ features extracted for MTFP ------------------%%%
    nucl_area_Nep                  =zeros(no_obj,numbM);
    cyt_area_Nep                   =zeros(no_obj,numbM);
    mean_Nep_int_per_area_C        =zeros(no_obj,numbM);
    mean_Nep_int_per_area_N        =zeros(no_obj,numbM);
    Nep_Conc_T                     =zeros(no_obj,numbM);
    Nep_Conc_C                     =zeros(no_obj,numbM);
    Nep_Conc_N                     =zeros(no_obj,numbM);
    Nep_mean_int_N                 =zeros(no_obj,numbM);
    Nep_mean_int_C                 =zeros(no_obj,numbM);
    tot_Nep_fl_cyt                 = zeros(no_obj,numbM);
    tot_Nep_fl_nuc                 = zeros(no_obj,numbM);  
        
   end
    
    %%%%%%%-------------- mCyOFP1 feature allocation ---------------%%%%%%%
   
    if fl_used(1,6)==1

            Cdc10_stds=zeros(no_obj,numbM);
            Cdc10_means=zeros(no_obj,numbM);
            Cdc10_status=zeros(1,no_obj);
    end
    
    
    tmp_sizes=zeros(1,no_obj); %allocation of segmented cells 
    Lcells1=all_obj.cells(:,:,numbM);
    parfor i5=1:no_obj
        tmp_sizes(i5)=tmp_sizes(i5)+sum(sum(Lcells1==i5));
    end
    max_allowed_cell_size=max_size_vs_largest_cell*max(tmp_sizes);
    s=round(sqrt(max_allowed_cell_size))+50; %?????????????????this will determine the size of the matrix we will put put_vac etc into
    
    for c_time=numbM:-1:1 % loop to go through all time points 
        current_time=c_time;  
        disp(c_time)
        Lcells=all_obj.cells(:,:,c_time); %a broadcast variable
        location_cell=zeros(no_obj,4);
             
 %%%---------- load images of the target channels        
        image_number=sprintf('%09d',c_time);
        
        if fl_used(1,1)==1
            IT=imread([path_h exp_name '/Pos' num2str(pos) '/img_' image_number '_470 mTFP_000.tif']);
            IT=double(IT);
            IT=medfilt2(IT,'symmetric');
            backgr_T=((double(IT)+1).*(~Lcells));
            backgr_T=sort(backgr_T(:));
            [tmpV,posH]=max(backgr_T>0);
            backgr_T=median(backgr_T(max([1 posH-1]):end))-1;
            IT=IT-backgr_T;
            all_obj.med_T(1,c_time)=backgr_T;
        end
        
        if fl_used(1,2)==1
            IG=imread([path_h exp_name '/Pos' num2str(pos) '/img_' image_number '_505 mNG_000.tif']);
            IG=double(IG);
            IG=medfilt2(IG,'symmetric');
            backgr_G=(((IG)+1).*(~Lcells));
            backgr_G=sort(backgr_G(:));
            [tmpV,posH]=max(backgr_G>0);
            backgr_G=median(backgr_G(max([1 posH-1]):end))-1;
            IG=IG-backgr_G;
            all_obj.med_G(1,c_time)  =backgr_G;
        end
        
        if fl_used(1,3)==1
            IK=imread([path_h exp_name '/Pos' num2str(pos) '/img_' image_number '_555 mKok_000.tif']);
            IK=double(IK);
            IK=medfilt2(IK,'symmetric');
            backgr_K=((double(IK)+1).*(~Lcells));
            backgr_K=sort(backgr_K(:));
            [tmpV,posH]=max(backgr_K>0);
            backgr_K=median(backgr_K(max([1 posH-1]):end))-1;
            IK=IK-backgr_K;
            all_obj.med_K(1,c_time)  =backgr_K;
        end
        
        if fl_used(1,4)==1
            IRu=imread([path_h exp_name '/Pos' num2str(pos) '/img_' image_number '_555 mRuby3_000.tif']);
            IRu=double(IRu);
            IRu=medfilt2(IRu,'symmetric');
            backgr_Ru=((double(IRu)+1).*(~Lcells));
            backgr_Ru=sort(backgr_Ru(:));
            [~,posH]=max(backgr_Ru>0);
            backgr_Ru=median(backgr_Ru(max([1 posH-1]):end))-1;
            IRu=(IRu-backgr_Ru);
            all_obj.med_Ru(1,c_time)  =backgr_Ru;
        end

        
        if fl_used(1,5)==1
            ICa=imread([path_h exp_name '/Pos' num2str(pos) '/img_' image_number '_615 nm_000.tif']);
            ICa=double(ICa);
            ICa=medfilt2(ICa,'symmetric');
            backgr_mCa=((double(ICa)+1).*(~Lcells));
            backgr_mCa=sort(backgr_mCa(:));
            [~,posH]=max(backgr_mCa>0);
            backgr_mCa=median(backgr_mCa(max([1 posH-1]):end))-1;
            ICa=(ICa-backgr_mCa);
            all_obj.med_mCa(1,c_time)  =backgr_mCa;
        end
        
       if fl_used(1,6)==1
           channel= '470 GFP_000';
           I = get_image_name_OAM_1(Imagepath,current_time,channel,type);
           I = 2*(I-mean(I(:)));
           I = medfilt2(I,'symmetric');
       end
        
       
       
        parfor cell_no=1:no_obj % parallel loop to go through each cell
            if c_time>cell_exists(cell_no,2) % check that cells existed  at this time point
                
                ccell=(Lcells==cell_no); %a temporary variable
                appr_vol(cell_no,c_time)  = appr_vol(cell_no,c_time)+(sum(ccell(:))).^(1.5); % cell volume is approximated as the 2D area exponentiated to (3/2) 
                
                %Get Cell coordinates in the  image
                cell_margin=1;
                [x_cn,y_cn]=get_cell_coord(ccell,cell_margin);
                location_cell(cell_no,:)=location_cell(cell_no,:)+[y_cn(1) y_cn(end) x_cn(1) x_cn(end)];

                
%%%%-------------------- feature extraction for mTFP ------------------%%%%
    
                 if fl_used(1,1)==1
                    
                    put_IT=double(Lcells==cell_no).*IT;
                    max_nucl_int_IT(cell_no,c_time)=max_nucl_int_IT(cell_no,c_time)+max(put_IT(:));
                    mean_IT_int_per_area_T(cell_no,c_time)=mean_IT_int_per_area_T(cell_no,c_time)+sum(double(put_IT(:)))./sum(double(put_IT(:)>0));
                    IT_Conc_T(cell_no,c_time) =IT_Conc_T(cell_no,c_time)+sum(put_IT(put_IT(:)>0))./appr_vol(cell_no,c_time);
                    [p_nuc_IT] =Gaussian_nuclear_fit_OAM(IT,put_IT,peak_cutoff,x_size,y_size,Lcells,cell_no);
                    p_cyt_IT =bwmorph((Lcells==cell_no).*bwmorph(~p_nuc_IT,'thicken',1),'erode',2);
                    nuc_area_IT_tmp=sum(p_nuc_IT(:));
                    cyt_area_MTFP_tmp=sum(p_cyt_IT(:));
                    nucl_area_TFP(cell_no,c_time)  = nucl_area_TFP(cell_no,c_time) + nuc_area_IT_tmp;
                    cyt_area_TFP(cell_no,c_time)   = cyt_area_TFP(cell_no,c_time) +cyt_area_MTFP_tmp;
                    tot_IT_fl_cyt = sum(sum(double(p_cyt_IT).*double(IT)));
                    tot_TFP_fl_cyt(cell_no,c_time)=tot_TFP_fl_cyt(cell_no,c_time)+ tot_IT_fl_cyt;
                    tot_IT_fl_nuc = sum(sum(double(p_nuc_IT).*double(IT)));
                    tot_TFP_fl_nuc(cell_no,c_time)=tot_TFP_fl_nuc(cell_no,c_time)+ tot_IT_fl_nuc;
                    mean_TFP_int_per_area_C(cell_no,c_time)=mean_TFP_int_per_area_C(cell_no,c_time)+tot_IT_fl_cyt./double(cyt_area_MTFP_tmp);
                    mean_TFP_int_per_area_N(cell_no,c_time)=mean_TFP_int_per_area_N(cell_no,c_time)+tot_IT_fl_nuc./double(nuc_area_IT_tmp);
                    nuc_V_IT=nuc_area_IT_tmp.^(1.5);
                    cyt_V_IT=appr_vol(cell_no,c_time)-nuc_area_IT_tmp.^(1.5);
                    TFP_Conc_T(cell_no,c_time)     =TFP_Conc_T(cell_no,c_time)+sum(double(put_IT(:)))./appr_vol(cell_no,c_time);
                    TFP_Conc_C(cell_no,c_time)     =TFP_Conc_C(cell_no,c_time)+tot_IT_fl_cyt./cyt_V_IT;
                    TFP_Conc_N(cell_no,c_time)     =TFP_Conc_N(cell_no,c_time) +tot_IT_fl_nuc./nuc_V_IT;
                    ccell_MTFP=double(ccell(y_cn,x_cn)).*double(IT(y_cn,x_cn));
                    put_MTFP=(ccell_MTFP>(IT_mean_modifier.*mean(ccell_MTFP(ccell_MTFP>0))));
                    put_MTFP=bwareaopen(put_MTFP,5,4);
                    put_MTFP=imfill(put_MTFP,'holes');
                    TFP_mean_int_N(cell_no,c_time)   = TFP_mean_int_N(cell_no,c_time) +sum(sum(put_MTFP.*ccell_MTFP))./sum(put_MTFP(:));
                    no_MTFP=ccell_MTFP.*double(~put_MTFP);
                    TFP_mean_int_C(cell_no,c_time) = TFP_mean_int_C(cell_no,c_time)+sum(no_MTFP(no_MTFP>0))./sum(no_MTFP(no_MTFP>0)>0);
                    M=false(s,s);
                    M(1:length(y_cn),1:length(x_cn))=logical(put_MTFP);

                 end %  fl_used
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Neon Green
                
                            
                 if fl_used(1,2)==1
                    
                    put_Green=double(Lcells==cell_no).*IG;
                    max_nucl_int_Green(cell_no,c_time)=max_nucl_int_Green(cell_no,c_time)+max(put_Green(:));
                    mean_Green_int_per_area_T(cell_no,c_time)=mean_Green_int_per_area_T(cell_no,c_time)+sum(double(put_Green(:)))./sum(double(put_Green(:)>0));
                    Green_Conc_T(cell_no,c_time) =Green_Conc_T(cell_no,c_time)+sum(put_Green(put_Green(:)>0))./appr_vol(cell_no,c_time);
                    [p_nuc_Green ] =Gaussian_nuclear_fit_OAM(IG,put_Green,peak_cutoff,x_size,y_size,Lcells,cell_no);
                    p_cyt_Green  =bwmorph((Lcells==cell_no).*bwmorph(~p_nuc_Green ,'thicken',1),'erode',2);
                    nuc_area_Green_tmp=sum(p_nuc_Green (:));
                    cyt_area_Neon_Green_tmp=sum(p_cyt_Green (:));
                   nucl_area_Neon_Green (cell_no,c_time)  = nucl_area_Neon_Green (cell_no,c_time) + nuc_area_Green_tmp;
                    cyt_area_Neon_Green (cell_no,c_time)   = cyt_area_Neon_Green (cell_no,c_time) +cyt_area_Neon_Green_tmp;
                    tot_Green_fl_cyt = sum(sum(double(p_cyt_Green ).*double(IG)));
                    tot_Neon_Green_fl_cyt(cell_no,c_time)=tot_Neon_Green_fl_cyt(cell_no,c_time)+ tot_Green_fl_cyt;
                    tot_Green_fl_nuc = sum(sum(double(p_nuc_Green ).*double(IG)));
                    tot_Neon_Green_fl_nuc(cell_no,c_time)=tot_Neon_Green_fl_nuc(cell_no,c_time)+ tot_Green_fl_nuc;
                    mean_Neon_Green_int_per_area_C(cell_no,c_time)=mean_Neon_Green_int_per_area_C(cell_no,c_time)+tot_Green_fl_cyt./double(cyt_area_Neon_Green_tmp);
                    mean_Neon_Green_int_per_area_N(cell_no,c_time)=mean_Neon_Green_int_per_area_N(cell_no,c_time)+tot_Green_fl_nuc./double(nuc_area_Green_tmp);
                    nuc_V_Green =nuc_area_Green_tmp.^(1.5);
                    cyt_V_Green =appr_vol(cell_no,c_time)-nuc_area_Green_tmp.^(1.5);
                    Neon_Green_Conc_T(cell_no,c_time)     =Neon_Green_Conc_T(cell_no,c_time)+sum(double(put_Green(:)))./appr_vol(cell_no,c_time);
                    Neon_Green_Conc_C(cell_no,c_time)     =Neon_Green_Conc_C(cell_no,c_time)+tot_Green_fl_cyt./cyt_V_Green ;
                    Neon_Green_Conc_N(cell_no,c_time)     =Neon_Green_Conc_N(cell_no,c_time) +tot_Green_fl_nuc./nuc_V_Green ;
                    ccell_Neon_Green =double(ccell(y_cn,x_cn)).*double(IG(y_cn,x_cn));
                    put_Neon_Green =(ccell_Neon_Green >(G_mean_modifier.*mean(ccell_Neon_Green (ccell_Neon_Green >0))));
                    put_Neon_Green =bwareaopen(put_Neon_Green ,5,4);
                    put_Neon_Green =imfill(put_Neon_Green ,'holes');
                    Neon_Green_mean_int_N(cell_no,c_time)   = Neon_Green_mean_int_N(cell_no,c_time) +sum(sum(put_Neon_Green .*ccell_Neon_Green ))./sum(put_Neon_Green (:));
                    no_Neon_Green =ccell_Neon_Green .*double(~put_Neon_Green );
                    Neon_Green_mean_int_C(cell_no,c_time) = Neon_Green_mean_int_C(cell_no,c_time)+sum(no_Neon_Green (no_Neon_Green >0))./sum(no_Neon_Green (no_Neon_Green >0)>0);
                    M=false(s,s);
                    M(1:length(y_cn),1:length(x_cn))=logical(put_Neon_Green );

                 end %  fl_used

                             
                
                
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%KOK
                if fl_used(1,3)==1
                    
                    put_K=double(Lcells==cell_no).*IK;
                    max_nucl_int_K(cell_no,c_time)=max_nucl_int_K(cell_no,c_time)+max(put_K(:));
                    mean_K_int_per_area_T(cell_no,c_time)=mean_K_int_per_area_T(cell_no,c_time)+sum(double(put_K(:)))./sum(double(put_K(:)>0));
                    K_Conc_T(cell_no,c_time) =K_Conc_T(cell_no,c_time)+sum(put_K(put_K(:)>0))./appr_vol(cell_no,c_time);
                    [p_nuc_K] =Gaussian_nuclear_fit_OAM(IK,put_K,peak_cutoff,x_size,y_size,Lcells,cell_no);
                    p_cyt_K =bwmorph((Lcells==cell_no).*bwmorph(~p_nuc_K,'thicken',1),'erode',2);
                    nuc_area_K_tmp=sum(p_nuc_K(:));
                    cyt_area_KOK_tmp=sum(p_cyt_K(:));
                    nucl_area_KOK(cell_no,c_time)  = nucl_area_KOK(cell_no,c_time) + nuc_area_K_tmp;
                    cyt_area_KOK(cell_no,c_time)   = cyt_area_KOK(cell_no,c_time) +cyt_area_KOK_tmp;
                    tot_K_fl_cyt = sum(sum(double(p_cyt_K).*double(IK)));
                    tot_KOK_fl_cyt(cell_no,c_time)=tot_KOK_fl_cyt(cell_no,c_time)+ tot_K_fl_cyt;
                    tot_K_fl_nuc = sum(sum(double(p_nuc_K).*double(IK)));
                    tot_KOK_fl_nuc(cell_no,c_time)=tot_KOK_fl_nuc(cell_no,c_time)+ tot_K_fl_nuc;
                    mean_KOK_int_per_area_C(cell_no,c_time)=mean_KOK_int_per_area_C(cell_no,c_time)+tot_K_fl_cyt./double(cyt_area_KOK_tmp);
                    mean_KOK_int_per_area_N(cell_no,c_time)=mean_KOK_int_per_area_N(cell_no,c_time)+tot_K_fl_nuc./double(nuc_area_K_tmp);
                    nuc_V_K=nuc_area_K_tmp.^(1.5);
                    cyt_V_K=appr_vol(cell_no,c_time)-nuc_area_K_tmp.^(1.5);
                    KOK_Conc_T(cell_no,c_time)     =KOK_Conc_T(cell_no,c_time)+sum(double(put_K(:)))./appr_vol(cell_no,c_time);
                    KOK_Conc_C(cell_no,c_time)     =KOK_Conc_C(cell_no,c_time)+tot_K_fl_cyt./cyt_V_K;
                    KOK_Conc_N(cell_no,c_time)     =KOK_Conc_N(cell_no,c_time) +tot_K_fl_nuc./nuc_V_K;
                    ccell_KOK=double(ccell(y_cn,x_cn)).*double(IK(y_cn,x_cn));
                    put_KOK=(ccell_KOK>(K_mean_modifier.*mean(ccell_KOK(ccell_KOK>0))));
                    put_KOK=bwareaopen(put_KOK,5,4);
                    put_KOK=imfill(put_KOK,'holes');
                    KOK_mean_int_N(cell_no,c_time)   = KOK_mean_int_N(cell_no,c_time) +sum(sum(put_KOK.*ccell_KOK))./sum(put_KOK(:));
                    no_KOK=ccell_KOK.*double(~put_KOK);
                    KOK_mean_int_C(cell_no,c_time) = KOK_mean_int_C(cell_no,c_time)+sum(no_KOK(no_KOK>0))./sum(no_KOK(no_KOK>0)>0);
                    M=false(s,s);
                    M(1:length(y_cn),1:length(x_cn))=logical(put_KOK);

                end %  fl_used
    
    
               if fl_used(1,4)==1
                    
                    put_Ru=double(Lcells==cell_no).*IRu;
                    max_nucl_int_Ru(cell_no,c_time)=max_nucl_int_Ru(cell_no,c_time)+max(put_Ru(:));
                    mean_Ru_int_per_area_T(cell_no,c_time)=mean_Ru_int_per_area_T(cell_no,c_time)+sum(double(put_Ru(:)))./sum(double(put_Ru(:)>0));
                    Ru_Conc_T(cell_no,c_time) =Ru_Conc_T(cell_no,c_time)+sum(put_Ru(put_Ru(:)>0))./appr_vol(cell_no,c_time);
                    [p_nuc_Ru] =Gaussian_nuclear_fit_OAM(IRu,put_Ru,peak_cutoff,x_size,y_size,Lcells,cell_no);
                    p_cyt_Ru =bwmorph((Lcells==cell_no).*bwmorph(~p_nuc_Ru,'thicken',1),'erode',2);
                    nuc_area_Ru_tmp=sum(p_nuc_Ru(:));
                    cyt_area_R555R_tmp=sum(p_cyt_Ru(:));
                     nucl_area_Ruby(cell_no,c_time)  = nucl_area_Ruby(cell_no,c_time) + nuc_area_Ru_tmp;
                    cyt_area_Ruby(cell_no,c_time)   = cyt_area_Ruby(cell_no,c_time) +cyt_area_R555R_tmp;
                    tot_Ru_fl_cyt = sum(sum(double(p_cyt_Ru).*double(IRu)));
                    tot_Ruby_fl_cyt(cell_no,c_time)=tot_Ruby_fl_cyt(cell_no,c_time)+ tot_Ru_fl_cyt;
                    tot_Ru_fl_nuc = sum(sum(double(p_nuc_Ru).*double(IRu)));
                    tot_Ruby_fl_nuc(cell_no,c_time)=tot_Ruby_fl_nuc(cell_no,c_time)+ tot_Ru_fl_nuc;
                    mean_Ruby_int_per_area_C(cell_no,c_time)=mean_Ruby_int_per_area_C(cell_no,c_time)+tot_Ru_fl_cyt./double(cyt_area_R555R_tmp);
                    mean_Ruby_int_per_area_N(cell_no,c_time)=mean_Ruby_int_per_area_N(cell_no,c_time)+tot_Ru_fl_nuc./double(nuc_area_Ru_tmp);
                    nuc_V_Ru=nuc_area_Ru_tmp.^(1.5);
                    cyt_V_Ru=appr_vol(cell_no,c_time)-nuc_area_Ru_tmp.^(1.5);
                    Ruby_Conc_T(cell_no,c_time)     =Ruby_Conc_T(cell_no,c_time)+sum(double(put_Ru(:)))./appr_vol(cell_no,c_time);
                    Ruby_Conc_C(cell_no,c_time)     =Ruby_Conc_C(cell_no,c_time)+tot_Ru_fl_cyt./cyt_V_Ru;
                    Ruby_Conc_N(cell_no,c_time)     =Ruby_Conc_N(cell_no,c_time) +tot_Ru_fl_nuc./nuc_V_Ru;
                    ccell_R555R=double(ccell(y_cn,x_cn)).*double(IRu(y_cn,x_cn));
                    put_R555R=(ccell_R555R>(R_mean_modifier.*mean(ccell_R555R(ccell_R555R>0))));
                    put_R555R=bwareaopen(put_R555R,5,4);
                    put_R555R=imfill(put_R555R,'holes');
                    Ruby_mean_int_N(cell_no,c_time)   = Ruby_mean_int_N(cell_no,c_time) +sum(sum(put_R555R.*ccell_R555R))./sum(put_R555R(:));
                    no_R555R=ccell_R555R.*double(~put_R555R);
                    Ruby_mean_int_C(cell_no,c_time) = Ruby_mean_int_C(cell_no,c_time)+sum(no_R555R(no_R555R>0))./sum(no_R555R(no_R555R>0)>0);
                    M=false(s,s);
                    M(1:length(y_cn),1:length(x_cn))=logical(put_R555R);

               end %  fl_used
        
        
                
                
                 if fl_used(1,5)==1
                    
                    put_Nep=double(Lcells==cell_no).*ICa;
                    max_nucl_int_Nep(cell_no,c_time)=max_nucl_int_Nep(cell_no,c_time)+max(put_Nep(:));
                    mean_Nep_int_per_area_T(cell_no,c_time)=mean_Nep_int_per_area_T(cell_no,c_time)+sum(double(put_Nep(:)))./sum(double(put_Nep(:)>0));
                    Nept_Conc_T(cell_no,c_time) =Nept_Conc_T(cell_no,c_time)+sum(put_Nep(put_Nep(:)>0))./appr_vol(cell_no,c_time);
                    [p_nuc_Nep] =Gaussian_nuclear_fit_OAM(ICa,put_Nep,peak_cutoff,x_size,y_size,Lcells,cell_no);
                    p_cyt_Nep =bwmorph((Lcells==cell_no).*bwmorph(~p_nuc_Nep,'thicken',1),'erode',2);
                    nuc_area_Nep_tmp=sum(p_nuc_Nep(:));
                    cyt_area_Nep_tmp=sum(p_cyt_Nep(:));
                    nucl_area_Nep(cell_no,c_time)  = nucl_area_Nep(cell_no,c_time) + nuc_area_Nep_tmp;
                    cyt_area_Nep(cell_no,c_time)   = cyt_area_Nep(cell_no,c_time) +cyt_area_Nep_tmp;
                    tot_N615N_fl_cyt = sum(sum(double(p_cyt_Nep).*double(ICa)));
                    tot_Nep_fl_cyt(cell_no,c_time)=tot_Nep_fl_cyt(cell_no,c_time)+ tot_N615N_fl_cyt;
                    tot_N615N_fl_nuc = sum(sum(double(p_nuc_Nep).*double(ICa)));
                    tot_Nep_fl_nuc(cell_no,c_time)=tot_Nep_fl_nuc(cell_no,c_time)+ tot_N615N_fl_nuc;
                    mean_Nep_int_per_area_C(cell_no,c_time)=mean_Nep_int_per_area_C(cell_no,c_time)+tot_N615N_fl_cyt./double(cyt_area_Nep_tmp);
                    mean_Nep_int_per_area_N(cell_no,c_time)=mean_Nep_int_per_area_N(cell_no,c_time)+tot_N615N_fl_nuc./double(nuc_area_Nep_tmp);
                    nuc_V_Nep=nuc_area_Nep_tmp.^(1.5);
                    cyt_V_Nep=appr_vol(cell_no,c_time)-nuc_area_Nep_tmp.^(1.5);
                    Nep_Conc_T(cell_no,c_time)     =Nep_Conc_T(cell_no,c_time)+sum(double(put_Nep(:)))./appr_vol(cell_no,c_time);
                    Nep_Conc_C(cell_no,c_time)     =Nep_Conc_C(cell_no,c_time)+tot_N615N_fl_cyt./cyt_V_Nep;
                    Nep_Conc_N(cell_no,c_time)     =Nep_Conc_N(cell_no,c_time) +tot_N615N_fl_nuc./nuc_V_Nep;
                    ccell_N615N=double(ccell(y_cn,x_cn)).*double(ICa(y_cn,x_cn));
                    put_N615N=(ccell_N615N>(FR_mean_modifier.*mean(ccell_N615N(ccell_N615N>0))));
                    put_N615N=bwareaopen(put_N615N,5,4);
                    put_N615N=imfill(put_N615N,'holes');
                    Nep_mean_int_N(cell_no,c_time)   = Nep_mean_int_N(cell_no,c_time) +sum(sum(put_N615N.*ccell_N615N))./sum(put_N615N(:));
                    no_N615N=ccell_N615N.*double(~put_N615N);
                    Nep_mean_int_C(cell_no,c_time) = Nep_mean_int_C(cell_no,c_time)+sum(no_N615N(no_N615N>0))./sum(no_N615N(no_N615N>0)>0);
                    M=false(s,s);
                    M(1:length(y_cn),1:length(x_cn))=logical(put_N615N);

                 end %  fl_used
                
%% this it quantifycation of the standard deviation of Cdc10 at cell periphery as cell cycle marker 
                 if fl_used(1,6)==1
                     
                    Lcells1=all_obj.cells(:,:,c_time);
                    Cells01=find(Lcells1>=1);
                    Lcells1(Cells01)=1;
                    Lcells1= imdilate(Lcells1,[1 1 1; 1 1 1; 1 1 1]); % increases cell contour by one pixel
                    I2=I.*uint16(Lcells1);                     
                    cell_margin=2;
                    [x_cn,y_cn]=get_cell_coord(ccell,cell_margin);
                    Mem1=bwmorph(ccell(y_cn,x_cn),'remove'); % creates a 1 pixel belt base on cell contour
                    ccell4=double(I2(y_cn,x_cn)); 
                    c5=double(Mem1.*ccell4); % obtain the contour pixels of the cell
                   
                    pixels=(find(c5(:)>0)); % counts all pixels on the cell contour
                    Cdc10_stds(cell_no,c_time)=std(c5(pixels)) ; % calculates STD
                    Cdc10_means(cell_no,c_time)=mean(c5(pixels));% calculates mean
     
                 end %  fl_used   
                
            end %end if cell exists
            
        end %parfor
        
    end %time-loop
    
   %%%%%------------- alllocation of general features ----------------%%%%%
    
    all_obj.appr_vol =appr_vol;%
    
    %%%------ results for mTFP1 ------%%%
    if fl_used(1,1)==1
        all_obj.mean_T_int_per_area_T        =mean_IT_int_per_area_T;
        all_obj.T_Conc_T                     =IT_Conc_T;
        all_obj.max_nucl_int_IT              =max_nucl_int_IT;
    end
    %%%------ results for mNeonGreen ------%%%
    if fl_used(1,2)==1
        all_obj.mean_Green_int_per_area_T    =mean_Green_int_per_area_T;
        all_obj.Green_Conc_T                 =Green_Conc_T;
         all_obj.max_nucl_int_Green          =max_nucl_int_Green;   
    end
    %----------------------------------------------------------------------    
   
    %%%------ results for mKOk ------%%%
    if fl_used(1,3)==1
        all_obj.mean_K_int_per_area_T        =mean_K_int_per_area_T;
        all_obj.K_Conc_T                     =K_Conc_T;
        all_obj.max_nucl_int_K               =max_nucl_int_K;
    end
     %---------------------------------------------------------------------
     
     %%%------ results for mRuby3 or Scarlet-I ------%%%
    if fl_used(1,4)==1
                      
        all_obj.mean_Ru_int_per_area_T      =mean_Ru_int_per_area_T;        
        all_obj.Ru_Conc_T                   =Ru_Conc_T;                   
        all_obj.max_nucl_int_Ru             =max_nucl_int_Ru;          
        
    end
    %---------------------------------------------------------------------
    
    
     %%%------ results for mNeptune2.5 ------%%%
    if fl_used(1,5)==1
                       
        all_obj.mean_Nep_int_per_area_T     =mean_Nep_int_per_area_T;
        all_obj.Nept_Conc_T                 =Nept_Conc_T;
        all_obj.max_nucl_int_Nep            =max_nucl_int_Nep;
        
    end
   %-----------------------------------------------------------------------
   
     %%%------ results for Cdc10-mCyOFP1 ------%%%
    if fl_used(1,6)==1
         all_obj.Cdc10_STD     =Cdc10_stds;
         all_obj.Cdc10_Means   =Cdc10_means;
    end

   
   
     %%%------ specific features for mTFP1 ------%%%
     
    if fl_used(1,1)==1
    all_obj.nucl_area_MTFP                  =nucl_area_TFP;
    all_obj.cyt_area_MTFP                   =cyt_area_TFP;
    all_obj.mean_MTFP_int_per_area_C        =mean_TFP_int_per_area_C ;
    all_obj.mean_MTFP_int_per_area_N        =mean_TFP_int_per_area_N;
    all_obj.MTFP_Conc_T                     =TFP_Conc_T;
    all_obj.MTFP_Conc_C                     =TFP_Conc_C;
    all_obj.MTFP_Conc_N                     =TFP_Conc_N;
    all_obj.MTFP_mean_int_N                 =TFP_mean_int_N   ;
    all_obj.MTFP_mean_int_C                 =TFP_mean_int_C ;
    all_obj.tot_MTFP_fl_cyt                 =tot_TFP_fl_cyt;
    all_obj.tot_MTFP_fl_nuc                 =tot_TFP_fl_nuc;
    end
     
     %%%------ results for mNeonGreen ------%%%
     
    if fl_used(1,2)==1
    all_obj.nucl_area_Neon_Green                  =nucl_area_Neon_Green;
    all_obj.cyt_area_Neon_Green                   =cyt_area_Neon_Green;
    all_obj.mean_Neon_Green_int_per_area_C        =mean_Neon_Green_int_per_area_C ;
    all_obj.mean_Neon_Green_int_per_area_N        =mean_Neon_Green_int_per_area_N;
    all_obj.Neon_Green_Conc_T                     =Neon_Green_Conc_T;
    all_obj.Neon_Green_Conc_C                     =Neon_Green_Conc_C;
    all_obj.Neon_Green_Conc_N                     =Neon_Green_Conc_N;
    all_obj.Neon_Green_mean_int_N                 =Neon_Green_mean_int_N   ;
    all_obj.Neon_Green_mean_int_C                 =Neon_Green_mean_int_C ;
    all_obj.tot_Neon_Green_fl_cyt                 =tot_Neon_Green_fl_cyt;
    all_obj.tot_Neon_Green_fl_nuc                 =tot_Neon_Green_fl_nuc;
    end
     
     
     %%%------ results for mKOk ------%%%
     
    if fl_used(1,3)==1
    all_obj.nucl_area_KOK                  =nucl_area_KOK;
    all_obj.cyt_area_KOK                   =cyt_area_KOK;
    all_obj.mean_KOK_int_per_area_C        =mean_KOK_int_per_area_C ;
    all_obj.mean_KOK_int_per_area_N        =mean_KOK_int_per_area_N;
    all_obj.KOK_Conc_T                     =KOK_Conc_T;
    all_obj.KOK_Conc_C                     =KOK_Conc_C;
    all_obj.KOK_Conc_N                     =KOK_Conc_N;
    all_obj.KOK_mean_int_N                 =KOK_mean_int_N   ;
    all_obj.KOK_mean_int_C                 =KOK_mean_int_C ;
    all_obj.tot_KOK_fl_cyt                 =tot_KOK_fl_cyt;
    all_obj.tot_KOK_fl_nuc                 =tot_KOK_fl_nuc;
    end
    
     
    %%%------ results for mRuby3 or mScarlet-I  ------%%%
    if fl_used(1,4)==1
    all_obj.nucl_area_Ruby                  =nucl_area_Ruby;
    all_obj.cyt_area_Ruby                   =cyt_area_Ruby;
    all_obj.mean_Ruby_int_per_area_C        =mean_Ruby_int_per_area_C ;
    all_obj.mean_Ruby_int_per_area_N        =mean_Ruby_int_per_area_N;
    all_obj.Ruby_Conc_T                     =Ruby_Conc_T;
    all_obj.Ruby_Conc_C                     =Ruby_Conc_C;
    all_obj.Ruby_Conc_N                     =Ruby_Conc_N;
    all_obj.Ruby_mean_int_N                 =Ruby_mean_int_N   ;
    all_obj.Ruby_mean_int_C                 =Ruby_mean_int_C ;
    all_obj.tot_Ruby_fl_cyt                 =tot_Ruby_fl_cyt;
    all_obj.tot_Ruby_fl_nuc                 =tot_Ruby_fl_nuc;
    end
        
    %%%------ results for mNeptune2.5 ------%%%
    if fl_used(1,5)==1
    all_obj.nucl_area_Nep                  =nucl_area_Nep;
    all_obj.cyt_area_Nep                   =cyt_area_Nep;
    all_obj.mean_Nep_int_per_area_C        =mean_Nep_int_per_area_C ;
    all_obj.mean_Nep_int_per_area_N        =mean_Nep_int_per_area_N;
    all_obj.Nep_Conc_T                     =Nep_Conc_T;
    all_obj.Nep_Conc_C                     =Nep_Conc_C;
    all_obj.Nep_Conc_N                     =Nep_Conc_N;
    all_obj.Nep_mean_int_N                 =Nep_mean_int_N   ;
    all_obj.Nep_mean_int_C                 =Nep_mean_int_C ;
    all_obj.tot_Nep_fl_cyt                 =tot_Nep_fl_cyt;
    all_obj.tot_Nep_fl_nuc                 =tot_Nep_fl_nuc;
    end
    
    path_save=path_h;
    name1=[exp_name '/' exp_name '_pos' num2str(pos) '_extr_fl_data_test_translocator'];
    save(fullfile(path_save,name1), 'all_obj','peak_cutoff','numbM','-v7.3');
end
delete(gcp('nocreate'))


