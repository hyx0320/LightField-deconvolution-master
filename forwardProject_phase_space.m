function projection = forwardProject_phase_space(H,realspace)
%% =========================================================================
%% 函数功能：相空间前向投影 (Phase-space Forward Projection)
%% 作用：计算【单个特定视角】下，3D 样本在传感器上形成的 2D 投影图。
%%       它是相空间反卷积算法中，模拟光线传播的核心基础模块。
%% =========================================================================
%% 输入参数：
%  @H         : [3D 矩阵] 相空间下的点扩散函数 (PSF)。
%               维度依次为：(传感器高Y, 传感器宽X, 深度Z)
%               物理意义：注意！这里的 H 不再包含视角维度，因为它已经是
%               从原 5D 大矩阵中，被外层程序抽取出来的【单一视角 3D 光斑】。
%
%  @realspace : [3D 矩阵] 当前估计的三维样本体数据。
%               维度依次为：(高Y, 宽X, 深度Z)
%               物理意义：当前猜测的发光物体，数值代表发光强度。
%
%% 输出参数：
%  @projection : [2D 矩阵] 计算出的 2D 投影图。
%                维度依次为：(传感器高Y, 传感器宽X)
%                物理意义：模拟当前视角下，光线透过微透镜打在传感器上的光斑总和。
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

    % 1. 准备一张纯黑的 2D 画布，大小与样本的 XY 尺寸一致
    projection = zeros(size(realspace,1), size(realspace,2));
    
    % 2. 遍历三维样本的所有深度层 (Z轴)
    % 注意：这里没有视角XY遍历，因为在调用外层 (deconv_phase_space.m) 
    % 已经通过循环将 13x13 个视角拆开了，每次只传一个视角的数据进来。
    for z=1:size(realspace,3) 
        
        % 3. 物理光学晕染 (逐层卷积)
        % 抽出第 z 层的 2D 样本切片 realspace(:,:,z)
        % 抽出第 z 层的 2D 光斑切片 H(:,:,z)
        % 用 conv2 模拟这一层的光点散开成模糊光斑的物理过程
        % 'same' 保证卷积后的输出尺寸和原切片完全一致
        projection = projection + conv2(realspace(:,:,z), H(:,:,z), 'same');
    end
    
end