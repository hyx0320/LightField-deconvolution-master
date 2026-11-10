function Backprojection = backwardProject_LFM(Ht, projection)
%% =========================================================================
%% 函数功能：后向投影 (Backward Projection)
%% 作用：将 2D 平面的误差图，沿着光学逆路径反推回 3D 空间，得到 3D 修正量。
%% =========================================================================
%% 输入参数：
%  @Ht         : [5D 矩阵] 翻转后的点扩散函数 (Transposed PSF)。
%                维度依次为：(传感器高Y, 传感器宽X, 视角X, 视角Y, 深度Z)
%                物理意义：代表光线沿着原路"逆向"收拢的物理规律。
%
%  @projection : [2D 矩阵] 待反推的二维投影图 (通常是当前迭代的误差比率图)。
%                维度依次为：(传感器高Y, 传感器宽X)
%                物理意义：代表当前 2D 成像结果与真实照片的差异。
%
%% 输出参数：
%  @Backprojection : [3D 矩阵] 反推得到的 3D 修正系数矩阵。
%                    维度依次为：(传感器高Y, 传感器宽X, 深度Z)
%                    物理意义：告诉 3D 猜测体每个体素应该被放大或缩小的倍数。
%% =========================================================================
%% Robert's code
%% Reference:  Robert Prevede, Young-Gyu Yoon, Maximilian Hoffmann, Nikita Pak.etc. 
%% "Simultaneous whole-animal 3D imaging of neuronal activity using light-field microscopy " 
%% in Nature Methods VOL.11 NO.7|July 2014.
    
    % 1. 获取尺寸参数
    x3length = size(Ht,5); % 深度层数 (Z轴)
    Nnum = size(Ht,3);     % 视角分辨率 (微透镜阵列每个孔下方的像素数)
    
    % 2. 初始化 3D 输出矩阵 (用于存放反推回来的3D修正量)
    Backprojection = zeros(size(projection, 1), size(projection, 2), x3length);
    % 准备一张 2D 全零黑布，用于后续局部的临时计算
    zeroSlice = zeros(  size(projection,1) , size(projection, 2));
    
    % =========================================================================
    % 【第 1 层循环】遍历 3D 样本的所有深度层 (cc)
    % =========================================================================
    for cc=1:x3length,
        % 准备当前深度的 2D 累加画布 (存放当前层所有视角反推回来的总和)
        tempSliceBack = zeroSlice;
        
        % =====================================================================
        % 【第 2、3 层循环】遍历微透镜的所有视角 (aa, bb)
        % =====================================================================
        for aa=1:Nnum,
            for bb=1:Nnum,     
                
                % 3. 提取专属"逆推光斑" (注意这里是 Ht，代表翻转后的 PSF)
                % 物理意义：沿着光线散开的原路反推回去，使用翻转的印章。
                Hts = squeeze(Ht( :,: ,aa,bb,cc));
                
                % 准备当前视角的临时黑布
                tempSlice = zeroSlice;
                
                % 4. 微透镜视角采样 (与正向投影完全一致)
                % 用"筛子"从 2D 误差图中，抽取属于当前视角的误差像素。
                tempSlice( (aa:Nnum:end) , (bb:Nnum:end) ) = projection( (aa:Nnum:end) , (bb:Nnum:end) );
                
                % 5. 核心操作：反向物理晕染 (卷积)
                % 用逆推光斑 Hts 去卷 2D 误差，把平面的误差沿着光路"倒推"回当前深度的 3D 空间。
                % 把当前视角倒推出来的误差，累加到当前层的总修正画布上。
                tempSliceBack = tempSliceBack + conv2(tempSlice, Hts, 'same');   
            end
        end
        
        % 6. 把这一层所有视角汇总的 3D 修正量，存入最终输出的 3D 矩阵中
        Backprojection(:,:,cc) = Backprojection(:,:,cc) + tempSliceBack;
    end
end