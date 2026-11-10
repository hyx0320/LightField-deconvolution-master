function Output = deconv_phase_space_E(LF)
% 实验组 E：无螺旋更新 (标准行列扫描)
    maxIter = 20; 
    load('PSF/LF_psf.mat');
    Nnum = size(H,3);
    
    % 1. 保留 PSF 重排与倒序折叠
    IMGsize = size(H,1)-mod((size(H,1)-Nnum),2*Nnum);
    psf = zeros(IMGsize,IMGsize,Nnum,Nnum,size(H,5));
    for z=1:size(H,5) 
        LFtmp = zeros(IMGsize,IMGsize,Nnum,Nnum);
        index1 = round(size(H,1)/2)-fix(size(LFtmp,1)/2);
        index2 = round(size(H,1)/2)+fix(size(LFtmp,1)/2);
        for ii=1:size(H,3)
            for jj=1:size(H,4)
                LFtmp(:,:,ii,jj) = im_shift3(squeeze(H(index1:index2,index1:index2,ii,jj,z)),ii-((Nnum+1)/2), jj-(Nnum+1)/2);
            end
        end
        multiWDF = zeros(Nnum,Nnum,size(LFtmp,1)/size(H,3),size(LFtmp,2)/size(H,4),Nnum,Nnum);
        for i=1:size(H,3)
            for j=1:size(H,4)
                for a=1:size(LFtmp,1)/size(H,3)
                    for b=1:size(LFtmp,2)/size(H,4)
                        multiWDF(i,j,a,b,:,:) = squeeze( LFtmp( (a-1)*Nnum+i, (b-1)*Nnum+j, :, : ) );
                    end
                end
            end
        end
        WDF = zeros(size(LFtmp,1),size(LFtmp,2),Nnum,Nnum);
        for a=1:size(LFtmp,1)/size(H,3)
            for c=1:Nnum
                x = Nnum*a + 1 - c;
                for b=1:size(LFtmp,2)/size(H,4)
                    for d=1:Nnum
                        y = Nnum*b + 1 - d;
                        WDF(x,y,:,:) = squeeze(multiWDF(:,:,a,b,c,d));
                    end
                end
            end
        end
        psf(:,:,:,:,z) = WDF;   
    end
    psf_t = zeros(size(psf));
    for ii=1:Nnum
        for jj=1:Nnum
            for cc=1:size(psf,5)
                psf_t(:,:,ii,jj,cc) = rot90(squeeze(psf(:,:,ii,jj,cc)),2 );
            end
        end
    end   
    
    % 2. 光场图提取与插值
    LR_phase_space = zeros(size(LF,1)/Nnum,size(LF,2)/Nnum,Nnum,Nnum);
    for i=1:Nnum
        for j=1:Nnum
            LR_phase_space(:,:,i,j) = LF(i:Nnum:end,j:Nnum:end); 
        end
    end
    phase_space = zeros(size(LF,1),size(LF,2),Nnum,Nnum);
    for i=1:Nnum
        for j=1:Nnum
            phase_space(:,:,i,j) = imresize(LR_phase_space(:,:,i,j),[size(LF,1),size(LF,2)],'cubic');
        end
    end
    phase_space(phase_space<0) = 0;
    
    % 3. 低通滤波
    x=-fix(size(phase_space,1)/2):fix(size(phase_space,1)/2);
    [xx,yy]=meshgrid(x,x);
    cut_off_freq=10; 
    sub_aperture_filter=zeros(size(phase_space,1),size(phase_space,2));
    sub_aperture_filter(xx.^2+yy.^2 <= cut_off_freq^2)=1; 
    for i=1:Nnum
        for j=1:Nnum      
            phase_space_f_with_filter=fftshift(fft2(squeeze(phase_space(:,:,i,j)))).*sub_aperture_filter;
            phase_space_with_filter(:,:,i,j)=abs( ifft2(ifftshift(phase_space_f_with_filter)) );
        end
    end
    phase_space=phase_space_with_filter;
    phase_space(phase_space<0)=0;
    
    % 4. 计算权重
    weight=squeeze(sum(sum(sum(psf,1),2),5)); 
    weight(find(isnan(weight))) = 0;
    weight=weight./sum(weight(:));            
    weight=weight-min(weight(:));
    weight(weight<0.0035)=0;                  
    weight=weight.*80;                    
    
    % 5. 初始化
    Xguess=ones(size(LF,1),size(LF,2),size(psf,5));
    Xguess=Xguess./sum(Xguess(:)).*sum(LF(:));
    Htf=zeros(size(LF,1),size(LF,2),size(psf,5),Nnum,Nnum);
    for u=1:Nnum
        for v=1:Nnum
            if weight(u,v)==0
                continue;
            else
                Htf(:,:,:,u,v)= backwardProject_phase_space(squeeze(psf_t(:,:,u,v,:)), ones( size(LF) ) );
            end
        end
    end
    
    % 6. 【剔除】螺旋顺序，直接按行列扫描
    
    % 7. 主循环
    for i=1:maxIter
        tic;
        % !!! 核心改动：使用普通的 for 循环，而不是 index1/index2
        for u = 1:Nnum
            for v = 1:Nnum
                if weight(u,v)==0
                    continue;
                else
                    XguessBefore=Xguess;
                    HXguess=forwardProject_phase_space(squeeze(psf(:,:,u,v,:)), Xguess);            
                    errorEM=squeeze(phase_space(:,:,u,v))./HXguess;
                    errorEM(~isfinite(errorEM))=0;                
                    XguessCor = backwardProject_phase_space(squeeze(psf_t(:,:,u,v,:)),errorEM);               
                    Xguess_add=XguessBefore.*XguessCor./squeeze(Htf(:,:,:,u,v));
                    Xguess_add(find(isnan(Xguess_add)))=0;
                    Xguess_add(find(isinf(Xguess_add)))=0;
                    Xguess_add(Xguess_add<0)=0;
                    Xguess=Xguess_add.*weight(u,v)+(1-weight(u,v)).*XguessBefore;
                    Xguess(Xguess<0) = 0;              
                end
            end
        end
        ttime1=toc;
        disp(['E_NoSpiral iter ',num2str(i),' took ',num2str(ttime1),' secs']);
        
        % ===== 新增保存代码 =====
        temp_Xguess = Xguess; 
        save(['ReconResult_phase_space_deconv_E/Xguess_E_iter_', num2str(i), '.mat'], 'temp_Xguess', '-v7.3');
    end
    
    Output = Xguess(Nnum*2+1:end-Nnum*2,Nnum*2+1:end-Nnum*2,6:end-5);
    save('ReconResult_phase_space_deconv_E/Recon_E_NoSpiral.mat','Output');
end