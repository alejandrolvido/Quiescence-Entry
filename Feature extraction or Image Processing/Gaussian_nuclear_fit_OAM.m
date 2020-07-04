function [p_nuc]=Gaussian_nuclear_fit_EZGI(IG,put_XFP,peak_cutoff,x_size,y_size,Lcells,cell_no)


    p_nuc=zeros(x_size,y_size);
    %-----gaussian fit--------------------
    h = fspecial('average', 5);
    I_gf=  imfilter(IG,h,'replicate');
    Ignew=I_gf;
    nuc_tmp=Ignew.*double(Lcells==cell_no);
    
    Amp=max(max(nuc_tmp));
    s_f=regionprops(bwlabel(nuc_tmp==Amp),'Centroid'); %no need for bwlabel
    xtmpcenter=s_f(1).Centroid(1,1);
    ytmpcenter=s_f(1).Centroid(1,2);
    wind_size=10;
    xtmpcoordhere=max(1,xtmpcenter-wind_size):min(y_size,xtmpcenter+wind_size);
    ytmpcoordhere=max(1,ytmpcenter-wind_size):min(x_size,ytmpcenter+wind_size);
    [X, Y] = meshgrid(xtmpcoordhere,ytmpcoordhere);
    
    sX=size(X);
    sY=size(Y);
    ok_size=wind_size.*2+1;
    if sX(1,1)==ok_size && sX(1,2)==ok_size && sY(1,1)==ok_size && sY(1,2)==ok_size
        In_part=double(nuc_tmp(round(diag(Y)),round(diag(X))));
        edge_med=median([In_part(1,:) In_part(end,:) In_part(:,1)' In_part(:,end)']);
        x0=xtmpcenter;
        y0=ytmpcenter;
        mu=[x0,y0];
        best_fit_score=1e9;
        best_fit_no=[0 0];
        for jj=1:1:25
            for jj2=1:1:25
                limit_h=round(sqrt(jj.*jj2));
                for ii=-(limit_h-1):2:(limit_h-1)
                    Sigma = [jj -ii; -ii jj2];
                    F = mvnpdf([X(:) Y(:)],mu,Sigma);
                    F = reshape(F,length(X(:,1)),length(Y(:,1)));
                    Fmax=max(max(F));
                    Z=(( ((double(Amp)-edge_med)./Fmax)).*F+edge_med);
                    Z=Z.*(In_part>0);
                    tmp_score=sum(sum(abs((Z-In_part))));
                    if tmp_score<best_fit_score
                        best_fit_no=Sigma;%[ii jj];
                        best_fit_score=tmp_score;
                    end
                end
            end
        end
        %----------------- get the best solution------------------------
        %Sigma = [best_fit_no(1,2) -best_fit_no(1,1); -best_fit_no(1,1) best_fit_no(1,2)];
        %F = mvnpdf([X(:) Y(:)],mu,Sigma);
        F = mvnpdf([X(:) Y(:)],mu,best_fit_no);
        F = reshape(F,length(X(:,1)),length(Y(:,1)));
        Fmax=max(max(F));
        Z=(( ((double(Amp)-edge_med)./Fmax)).*F+edge_med);
        Z=Z.*(In_part>0);
        %---------------------------------------------------------------    
        p_nuc(round(diag(Y)),round(diag(X)))=p_nuc(round(diag(Y)),round(diag(X)))+(Z>(double(Amp).*peak_cutoff));
        %---------------------------------------------------------------------     
        p_nuc=double(p_nuc); %
    end %size