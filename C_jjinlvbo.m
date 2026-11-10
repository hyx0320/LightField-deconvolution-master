function Output = deconv_phase_space_C(LF)
% 实验组 C：只有滤波 (无重排/倒序)
    maxIter = 20; 
    load('PSF/LF_psf.mat');
    Nnum = size(H,3);
    
    % 1. 【剔除】PSF重排与倒序折叠
    % 直接使用未对齐、未折叠的原始 5D PSF (退化)
    psf = H; 
    psf_t = Ht; % 假设 .mat 文件里自带了 Ht (原版代码 load 后就有 Ht)
    
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
    
    % 3. 保留低通滤波
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
    
    % 4. 计算权重 (此时 psf 是原始的 H)
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
    
    % 6. 螺旋顺序 (原样)
    index1=[1,1,1,1,1,1,1,1,1,1,1,1,1,2,3,4,5,6,7,8,9,10,11,12,13,13,13,13,13,13,13,13,13,13,13,13,13,12,11,10,9,8,7,6,5,4,3,2,2,2,2,2,2,2,2,2,2,2,2,3,4,5,6,7,8,9,10,11,12,12,12,12,12,12,12,12,12,12,12,11,10,9,8,7,6,5,4,3,3,3,3,3,3,3,3,3,3,4,5,6,7,8,9,10,11,11,11,11,11,11,11,11,11,10,9,8,7,6,5,4,4,4,4,4,4,4,4,5,6,7,8,9,10,10,10,10,10,10,10,9,8,7,6,5,5,5,5,5,5,6,7,8,9,9,9,9,9,8,7,6,6,6,6,7,8,8,8,7,7];
    index2=[1,2,3,4,5,6,7,8,9,10,11,12,13,13,13,13,13,13,13,13,13,13,13,13,13,12,11,10,9,8,7,6,5,4,3,2,1,1,1,1,1,1,1,1,1,1,1,1,2,3,4,5,6,7,8,9,10,11,12,12,12,12,12,12,12,12,12,12,12,11,10,9,8,7,6,5,4,3,2,2,2,2,2,2,2,2,2,2,3,4,5,6,7,8,9,10,11,11,11,11,11,11,11,11,11,10,9,8,7,6,5,4,3,3,3,3,3,3,3,3,4,5,6,7,8,9,10,10,10,10,10,10,10,9,8,7,6,5,4,4,4,4,4,4,5,6,7,8,9,9,9,9,9,8,7,6,5,5,5,5,6,7,8,8,8,7,6,6,7];
    
    % 7. 主循环
    for i=1:maxIter
        tic;
        for u_2=1:13
            for v_2=1:13
                u=index1((u_2-1)*13+v_2);
                v=index2((u_2-1)*13+v_2);
                if weight(u,v)==0
                    continue;
                else
                    XguessBefore=Xguess;
                    % 直接用未对齐的原始 H 切片进行前向投影
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
        disp(['C_NoRearrange_WithFilter iter ',num2str(i),' took ',num2str(ttime1),' secs']);

        % ===== 新增保存代码 =====
        temp_Xguess = Xguess; 
        save(['ReconResult_phase_space_deconv_C/Xguess_C_iter_', num2str(i), '.mat'], 'temp_Xguess', '-v7.3');
    end
    
    Output = Xguess(Nnum*2+1:end-Nnum*2,Nnum*2+1:end-Nnum*2,6:end-5);
    
    % ⚠️ 注意修改保存路径
    save('ReconResult_phase_space_deconv_C/Recon_C_NoRearrange_WithFilter.mat','Output');
end