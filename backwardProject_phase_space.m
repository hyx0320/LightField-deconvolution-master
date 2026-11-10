function Backprojection = backwardProject_phase_space(Ht,projection)
%% =========================================================================
%% 函数功能：相空间后向投影 (Phase-space Backward Projection)
%% 作用：将【单个特定视角】的 2D 误差图，沿着光学逆路径反推回 3D 空间，
%%       计算出当前视角下的 3D 修正量矩阵。
%% =========================================================================
%% 输入参数：
%  @Ht         : [3D 矩阵] 相空间下翻转后的点扩散函数 (Transposed PSF)。
%                维度依次为：(传感器高Y, 传感器宽X, 深度Z)
%                物理意义：代表沿着光路逆向收拢的光学规律（单视角）。
%
%  @projection : [2D 矩阵] 待反推的二维投影图 (当前视角下的误差比率图)。
%                维度依次为：(传感器高Y, 传感器宽X)
%                物理意义：代表当前视角下，模拟成像与真实拍到的子孔径图之间的差异。
%
%% 输出参数：
%  @Backprojection : [3D 矩阵] 反推得到的当前视角 3D 修正系数矩阵。
%                    维度依次为：(传感器高Y, 传感器宽X, 深度Z)
%                    物理意义：告诉 3D 猜测体在当前视角下，每一层每个体素该被放大或缩小的倍数。
%% =========================================================================
% The Code is created based on the method described in the following paper 
%   [1]  ZHI LU, JIAMIN WU, HUI QIAO and YOU ZHOU.etc,
%        Phase-space deconvolution for light field microscopy. 
%   [2] Robert Prevede, Young-Gyu Yoon, Maximilian Hoffmann, Nikita Pak.etc. 
%       "Simultaneous whole-animal 3D imaging of neuronal activity using light-field microscopy "   
%       in Nature Methods VOL.11 NO.7|July 2014.
%
%    Author: ZHI LU (luz18@mails.tsinghua.edu.cn)
%    Date  : 05/21/2019

    % 1. 初始化 3D 输出画布
    % 高度和宽度等于传入的 2D 误差图，深度等于 PSF 的层数。
    Backprojection = zeros(size(projection,1), size(projection,2), size(Ht,3));
    
    % 2. 遍历三维空间的每一个深度层 (z)
    for z=1:size(Ht,3)
        
        % 3. 核心操作：逐层反向卷积
        % 拿 2D 的误差图 projection，和当前深度的翻转光斑 Ht(:,:,z) 进行卷积。
        % 物理意义：把平面上的 2D 误差，沿着逆向光路"分摊"回第 z 层 3D 空间。
        % （这里没有视角循环和筛子操作，因为传进来的 Ht 和 projection 已经是单视角的了）
        Backprojection(:,:,z) = conv2(projection, Ht(:,:,z), 'same');
    end
end