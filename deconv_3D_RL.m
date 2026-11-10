function Output=deconv_3D_RL(LF)
%% =========================================================================
%% 函数功能：3D Richardson-Lucy 反卷积算法 (主循环)
%% 作用：已知一张真实的 2D 光场照片(LF)和光学系统的点扩散函数(PSF)，
%%       通过不断"模拟拍照 -> 算误差 -> 修正猜测"的迭代，反推出真实的 3D 样本。
%% =========================================================================
%% 输入参数：
%  @LF         : [2D 矩阵] 光场相机拍到的真实二维原始照片。
%
%% 输出参数：
%  @Output     : [3D 矩阵] 最终重建出的三维样本（已切除了边缘伪影）。
%% =========================================================================
%% Reference:  Robert Prevede, Young-Gyu Yoon, Maximilian Hoffmann, Nikita Pak.etc. 
%% "Simultaneous whole-animal 3D imaging of neuronal activity using light-field microscopy " 
%% in Nature Methods VOL.11 NO.7|July 2014.
%% Downloaded from https://www.nature.com/articles/nmeth.2964
%% =========================================================================
    % 1. 加载预先计算好的 PSF 图鉴 (包含正向 H 和 逆向 Ht)
    load('PSF/LF_psf.mat');
    
    % 2. 获取系统参数
    Nnum=size(H,3);     % 视角分辨率 (微透镜下覆盖的像素数，比如 13)
    maxIter=20;         % 最大迭代次数 (通常 10~30 次即可收敛)
    
    % 3. 初始化 3D 猜测体 (捏一个初始泥块)
    % 创建一个与真实照片 XY 尺寸相同，与 PSF 深度层数相同的全 1 矩阵
    Xguess=ones(size(LF,1),size(LF,2),size(H,5));
    
    % 4. 能量归一化 (对齐总发光量)
    % 把猜测体的总亮度，缩放到和真实照片 LF 的总亮度一样。
    Xguess=Xguess./sum(Xguess(:)).*sum(LF(:));
    
    % 5. 生成"终极参考底片" (Target)
    % 把真实的 2D 照片，沿着逆向光路反推回 3D，作为每次比对误差的基准。
    Htf=backwardProject_LFM(Ht,LF); 
    
    % =========================================================================
    % 开始迭代"雕刻"
    % =========================================================================
    for i=1:maxIter,
        tic; % 开始计时
        
        % 步骤 A：正演模拟（虚拟拍照）
        % 把当前猜的 3D 泥块，压扁成一张 2D 照片，看看长什么样。
        HXguess = forwardProject_LFM(H,Xguess); 
        
        % 步骤 B：反演模拟（准备对比）
        % 把刚才模拟出来的 2D 照片，沿着逆向光路反推回 3D 空间。
        HXguessBack = backwardProject_LFM(Ht,HXguess); 
        
        % 步骤 C：计算误差比率
        % 拿"真实底片"除以"模拟底片"。
        % 结果 > 1 说明该位置偏暗，结果 < 1 说明该位置偏亮。
        errorBack = Htf ./ HXguessBack; 
        
        % 步骤 D：乘性更新（雕刻泥块）
        % 根据误差比率，把原猜测体该亮的地方调亮，该暗的地方压暗。
        Xguess = Xguess .* errorBack; 
        
        % 步骤 E：安全清洗（物理约束）
        % 消除 0 除以 0 产生的 NaN (非数)，以及除以极小值产生的 Inf (无穷大)。
        Xguess(find(isnan(Xguess))) = 0;
        Xguess(find(isinf(Xguess))) = 0;
        % 强制归零负数，因为光子数量不可能为负。
        Xguess(Xguess<0) = 0;
        
        ttime1=toc; % 停止计时
        
        % 打印本轮迭代信息，方便监控进度
        disp(['iter ',num2str(i),' | ',num2str(maxIter),', 3D RL deconvolution took ',num2str(ttime1),' secs']);
    
    
        % ===== 可选：存档代码 =====
        % 注意：此段非原作者代码，是后来为了观察收敛过程而加上的。
        % 每轮迭代都会产生一个 ~45MB 的文件，正式跑数据时建议注释掉！
        % 将当前迭代的未裁剪结果保存下来
        temp_Xguess = Xguess; 
        save(['ReconResult_3D_RL_deconv/Xguess_iter_', num2str(i), '.mat'], 'temp_Xguess', '-v7.3');
    end
    
    % =========================================================================
    % 后期处理：切边裁剪
    % =========================================================================
    % 因为卷积的边界效应，图像四周会出现一圈扭曲的伪影。
    % 切掉四周 (Nnum*2 = 26) 个像素，切掉首尾各 5 层深度，只保留最干净的核心区。
    Output=Xguess(Nnum*2+1:end-Nnum*2,Nnum*2+1:end-Nnum*2,6:end-5);
end