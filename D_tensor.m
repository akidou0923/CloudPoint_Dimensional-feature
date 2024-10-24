function [ S ] = D_tensor( NB_points,j )
%   此函数根据给定3D点以及K个最近领域点计算3D结构张量S

cord_sum = sum(NB_points);      %计算点云数据的坐标和(sum默认按列求和)

G_center = cord_sum./j;    %求邻域的坐标平均值
G=repmat(G_center,j,1);    %构造由j个重心坐标行向量组成的矩阵，j为邻域点个数
c_matrix =(NB_points-G)'*(NB_points-G);

% c_matrix = zeros(3,3); % 协方差阵3*3
% 
% 
% NB_points(j,:)-G_center;
% 
% for i=1:j  %逐点计算协方差阵
%     x1 = NB_points(i,:) - G_center;
%     c_matrix = c_matrix + (x1'*x1); 
% end

S = c_matrix./j;        %计算三维结构张量S
end

