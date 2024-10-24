function [ S ] = D_tensor( NB_points,j )
%   �˺������ݸ���3D���Լ�K�������������3D�ṹ����S

cord_sum = sum(NB_points);      %����������ݵ������(sumĬ�ϰ������)

G_center = cord_sum./j;    %�����������ƽ��ֵ
G=repmat(G_center,j,1);    %������j������������������ɵľ���jΪ��������
c_matrix =(NB_points-G)'*(NB_points-G);

% c_matrix = zeros(3,3); % Э������3*3
% 
% 
% NB_points(j,:)-G_center;
% 
% for i=1:j  %������Э������
%     x1 = NB_points(i,:) - G_center;
%     c_matrix = c_matrix + (x1'*x1); 
% end

S = c_matrix./j;        %������ά�ṹ����S
end

