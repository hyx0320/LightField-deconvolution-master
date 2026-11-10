function Output=deconv_phase_space(LF)
%% =========================================================================
%% 函数功能：相空间反卷积主程序 (Phase-space Deconvolution)
%% 作用：基于相空间重构与 Ptychographic (叠层成像) 策略，将混叠的 2D 光场图
%%       重建为高分辨率的三维样本体数据。
%% =========================================================================
%% 输入参数：
%  @LF         : [2D 矩阵] 光场相机拍到的真实二维原始照片。
%
%% 输出参数：
%  @Output     : [3D 矩阵] 最终重建出的三维样本（已切除了边缘伪影）。
%% =========================================================================
%% Reference:  ZHI LU, JIAMIN WU, HUI QIAO and YOU ZHOU .etc,
%%             Phase-space deconvolution for light field microscopy.
%% Author: ZHI LU (luz18@mails.tsinghua.edu.cn)
%% Date  : 05/21/2019
%% =========================================================================
    %% 最大迭代次数
    maxIter=2; % 注意：原作者这里写了2，但内部有169次子循环，实际计算量很大。
    
    %% ================== 第一阶段：预处理 (数据重排) ==================
    load('PSF/LF_psf.mat');
    Nnum=size(H,3);
    
    % 1. 将 5D 的 H 重排为相空间 PSF
    IMGsize=size(H,1)-mod((size(H,1)-Nnum),2*Nnum);
    psf =zeros(IMGsize,IMGsize,Nnum,Nnum,size(H,5));
    
    for z=1:size(H,5) % 遍历深度
        LFtmp=zeros(IMGsize,IMGsize,Nnum,Nnum);
        index1=round(size(H,1)/2)-fix(size(LFtmp,1)/2);
        index2=round(size(H,1)/2)+fix(size(LFtmp,1)/2);
        
        % 核心操作：把每个微透镜下的子孔径光斑平移对齐到中心
        for ii=1:size(H,3)
            for jj=1:size(H,4)
                LFtmp(:,:,ii,jj)=im_shift3(squeeze(H(index1:index2,index1:index2,ii,jj,z)),ii-((Nnum+1)/2), jj-(Nnum+1)/2);
            end
        end
        
        % 复杂的降维重组，生成相空间(Wigner分布)格式的 PSF (multiWDF)
        multiWDF=zeros(Nnum,Nnum,size(LFtmp,1)/size(H,3),size(LFtmp,2)/size(H,4),Nnum,Nnum);
        for i=1:size(H,3)
            for j=1:size(H,4)
                for a=1:size(LFtmp,1)/size(H,3)
                    for b=1:size(LFtmp,2)/size(H,4)
                        multiWDF(i,j,a,b,:,:)=squeeze( LFtmp(  (a-1)*Nnum+i,(b-1)*Nnum+j,:,:  )  );
                    end
                end
            end
        end
        
        % 重排回物理坐标
        WDF=zeros(  size(LFtmp,1),size(LFtmp,2),Nnum,Nnum  );
        for a=1:size(LFtmp,1)/size(H,3)%微透镜个数
            for c=1:Nnum%视角个数
                x=Nnum*a+1-c;
                for b=1:size(LFtmp,2)/size(H,4)%微透镜个数
                    for d=1:Nnum%视角个数
                        y=Nnum*b+1-d;
                        WDF(x,y,:,:)=squeeze(multiWDF(:,:,a,b,c,d));
                    end
                end
            end
        end
        psf(:,:,:,:,z)=WDF;   
    end
    
    % 2. 翻转 PSF 生成 psf_t (用于后向投影)
    psf_t=zeros(size(psf));
    for ii=1:Nnum
        for jj=1:Nnum
            for cc=1:size(psf,5)
                psf_t(:,:,ii,jj,cc)=rot90(squeeze(psf(:,:,ii,jj,cc)),2 ); % 旋转180度
            end
        end
    end   
    
    % 3. 将真实的 2D 光场图重排为相空间 (子孔径图)
    LR_phase_space=zeros(size(LF,1)/Nnum,size(LF,2)/Nnum,Nnum,Nnum);
    for i=1:Nnum
        for j=1:Nnum
            % 用筛子抽出每个视角的子孔径图
            LR_phase_space(:,:,i,j)=LF(i:Nnum:end,j:Nnum:end); 
        end
    end
    
    % 插值放大回原图尺寸
    phase_space = zeros( size(LF,1),size(LF,2),Nnum,Nnum );
    for i = 1:Nnum
        for j = 1:Nnum
            phase_space(:,:,i,j) = imresize(LR_phase_space(:,:,i,j),[size(LF,1),size(LF,2)],'cubic');
        end
    end
    phase_space(phase_space<0)=0;
    
    % 4. 频域低通滤波 (减少视角间的串扰)
    x=-fix(size(phase_space,1)/2):fix(size(phase_space,1)/2);
    [xx,yy]=meshgrid(x,x);
    cut_off_freq=10; % 截止频率
    sub_aperture_filter=zeros(size(phase_space,1),size(phase_space,2));
    sub_aperture_filter(xx.^2+yy.^2 <= cut_off_freq^2)=1; % 生成圆形掩膜
    
    for i=1:Nnum
        for j=1:Nnum      
            % 傅里叶变换 -> 乘掩膜 -> 逆变换
            phase_space_f_with_filter=fftshift(fft2(squeeze(phase_space(:,:,i,j)))).*sub_aperture_filter;
            phase_space_with_filter(:,:,i,j)=abs( ifft2(ifftshift(phase_space_f_with_filter)) );
        end
    end
    phase_space=phase_space_with_filter;
    phase_space(phase_space<0)=0;
    
    % 5. 计算每个视角的权重 (边缘视角光线弱，权重低)
    weight=squeeze(sum(sum(sum(psf,1),2),5)); % 对 XY 和 Z 积分，得到 13x13 的视角能量矩阵
    weight(find(isnan(weight))) = 0;
    weight=weight./sum(weight(:));            % 归一化
    weight=weight-min(weight(:));
    weight(weight<0.0035)=0;                  % 剔除太弱的视角
    weight=weight.*80;                        % 放大权重数值
    
    %% ================== 第二阶段：初始化 ==================
    
    % 初始化 3D 猜测体 (同 RL 算法)
    Xguess=ones(size(LF,1),size(LF,2),size(psf,5));
    Xguess=Xguess./sum(Xguess(:)).*sum(LF(:));
    
    % 预计算每个视角的 3D 参考底片 Htf
    Htf=zeros( size(LF,1),size(LF,2) , size(psf,5) , Nnum ,Nnum );
    for u=1:Nnum
        for v=1:Nnum
            if weight(u,v)==0
                continue; % 跳过无效视角
            else
                % 注意：这里传入的是全 1 矩阵 ones(size(LF))，代表理论的均匀背景
                Htf(:,:,:,u,v)= backwardProject_phase_space(squeeze(psf_t(:,:,u,v,:)), ones( size(LF) ) );
            end
        end
    end
    
    % 定义 Ptychographic 螺旋扫描顺序 (从中心向四周扩散，共 13x13=169 个位置)
    index1=[1,1,1,1,1,1,1,1,1,1,1,1,1,2,3,4,5,6,7,8,9,10,11,12,13,13,13,13,13,13,13,13,13,13,13,13,13,12,11,10,9,8,7,6,5,4,3,2,2,2,2,2,2,2,2,2,2,2,2,3,4,5,6,7,8,9,10,11,12,12,12,12,12,12,12,12,12,12,12,11,10,9,8,7,6,5,4,3,3,3,3,3,3,3,3,3,3,4,5,6,7,8,9,10,11,11,11,11,11,11,11,11,11,10,9,8,7,6,5,4,4,4,4,4,4,4,4,5,6,7,8,9,10,10,10,10,10,10,10,9,8,7,6,5,5,5,5,5,5,6,7,8,9,9,9,9,9,8,7,6,6,6,6,7,8,8,8,7,7];
    index2=[1,2,3,4,5,6,7,8,9,10,11,12,13,13,13,13,13,13,13,13,13,13,13,13,13,12,11,10,9,8,7,6,5,4,3,2,1,1,1,1,1,1,1,1,1,1,1,1,2,3,4,5,6,7,8,9,10,11,12,12,12,12,12,12,12,12,12,12,12,11,10,9,8,7,6,5,4,3,2,2,2,2,2,2,2,2,2,2,3,4,5,6,7,8,9,10,11,11,11,11,11,11,11,11,11,10,9,8,7,6,5,4,3,3,3,3,3,3,3,3,4,5,6,7,8,9,10,10,10,10,10,10,10,9,8,7,6,5,4,4,4,4,4,4,5,6,7,8,9,9,9,9,9,8,7,6,5,5,5,5,6,7,8,8,8,7,6,6,7];
    
    %% ================== 第三阶段：核心迭代 (Ptychographic) ==================
    for i=1:maxIter
        tic;
        % 按照螺旋顺序，逐个视角更新 3D 体
        for u_2=1:13
            for v_2=1:13
                u=index1((u_2-1)*13+v_2);
                v=index2((u_2-1)*13+v_2);
                
                if weight(u,v)==0
                    continue; % 跳过无效视角
                else
                    XguessBefore=Xguess; % 备份更新前的状态
                    
                    % 1. 正演模拟：算出当前 3D 样本在 (u,v) 视角下的模拟投影图
                    HXguess=forwardProject_phase_space(squeeze(psf(:,:,u,v,:)), Xguess);            
                    
                    % 2. 计算误差比率 (真实该视角图 / 模拟图)
                    errorEM=squeeze(phase_space(:,:,u,v))./HXguess;
                    errorEM(~isfinite(errorEM))=0;                
                    
                    % 3. 反演模拟：把 2D 误差反推回 3D 空间
                    XguessCor = backwardProject_phase_space(squeeze(psf_t(:,:,u,v,:)),errorEM) ;  
                    
                    % 4. 计算当前视角下的 3D 更新提议值
                    Xguess_add=XguessBefore.*XguessCor./squeeze(Htf(:,:,:,u,v));
                    Xguess_add(find(isnan(Xguess_add))) = 0;
                    Xguess_add(find(isinf(Xguess_add))) = 0;
                    Xguess_add(Xguess_add<0 ) = 0;
                    
                    % 5. 加权融合 (关键)：按视角权重混合新旧 3D 体
                    % 权重越大，越相信本次更新；权重小，则保留更多旧状态
                    Xguess=Xguess_add.*weight(u,v)+(1-weight(u,v)).*XguessBefore;
                    Xguess(Xguess<0) = 0;              
                end
            end
        end
        
        ttime1=toc;
        disp(['iter ',num2str(i),' | ',num2str(maxIter),', phase-space deconvolution took ',num2str(ttime1),' secs']);
    
        % ===== 可选：存档代码 (跑大图时建议注释掉) =====
        temp_Xguess = Xguess; 
        save(['ReconResult_phase_space_deconv/Xguess_iter_', num2str(i), '.mat'], 'temp_Xguess', '-v7.3');
    end
    
    %% ================== 第四阶段：后处理 ==================
    % 切除边缘伪影 (受卷积边界效应影响的部分)
    Output=Xguess(Nnum*2+1:end-Nnum*2,Nnum*2+1:end-Nnum*2,6:end-5);
end