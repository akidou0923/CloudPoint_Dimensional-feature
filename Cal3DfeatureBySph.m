clc,clear;
tic;

%%%%%%%%%%%%%%%%%%%此代码用来计算点云的三维特征(球形邻域），并给出最佳邻域半径（单位为m)%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%输入一个点云，输出该点云在最佳邻域半径下的planarity（平面性）、linearity（线性）以及sphericity（立体性）

Point_data = load ('test.txt'); %读取点云数据
MD1 = KDTreeSearcher(Point_data); %建立KD-tree模型
[row,line] = size(Point_data);%%（原始点的行列）

C_Lambda_all = zeros(row,1);
D_Feature_all = zeros(row,3);
% E_Dim_all = zeros(row,1);

rSize=zeros(10,1);
rSize(1)=0.1;

for i= 2:10
    rSize(i)=0.1+(i^2)*0.03;
end

toc

 for i = 1:row %遍历每点
     tic
    C_Lambda = zeros(10,1); 
    D_Feature = zeros(10,3);%用来存储该点不同邻域的LPS
    E_Dim = zeros(10,1);    %用来存储该点不同邻域的香农熵
  
    for k = 1:10 
        j=Size(k,1);
        Idx = rangesearch(MD1,Point_data(i,:),j);%从Point_data中根据邻域半径找临近点，返回索引
        NB_point = Point_data(Idx{1}',:);%从Point_data中根据索引得到临近点集

        [r,l]=size(NB_point);
        S = D_tensor(NB_point,r);  %%计算协方差阵
        [Z,L] = eig(S);          %%特征值
        Normal_v(1,:) = (Z(:,1))';   %将法向量存在矩阵Normal_v中

        [d,ind]=sort(diag(L),'descend'); %%将特征值按降序的方法排列
        Lambda =d';

% %%%%%%%论文里的方法，对特征值取根号%%%%%%%%%%%%%%
%          e1=sqrt(Lambda(1));
%          e2=sqrt(Lambda(2));
%          e3=sqrt(Lambda(3));
        
                
%%%%%%%%%%%%%%%%%%%% 原始代码，用来特征值归一化%%%%%%%%%%%%%%
        lambda_sum = sum(Lambda);
        e1 = Lambda(1)/lambda_sum;
        e2 = Lambda(2)/lambda_sum;
        e3 = Lambda(3)/lambda_sum;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        L_lambda = (e1-e2) / e1;
        P_lambda = (e2-e3) /e1;
        S_lambda = e3  /e1 ; %计算维度特征L_lambda、P_lambda、S_lambda
        Feature = [L_lambda,P_lambda,S_lambda]; 
        
        E_dim = -L_lambda*log(L_lambda)-P_lambda*log(P_lambda)-S_lambda*log(S_lambda);
        D_Feature(k,:) = Feature;
        E_Dim(k,1) = E_dim;
    end
        [~,OT] = min(E_Dim);   %找到最小的香农熵对应的下标
        D_Feature_all(i,:) = D_Feature(OT,:);
%       E_Dim_all(i,:) = E_Dim(OT,:);
%         
  if mod(i,1000) == 0 
        fprintf('第%d个点\n',i);
 end
toc
 end
 
All_feature = [Point_data,D_Feature_all];
% fid = fopen('feature_3D.txt','wt');
% for n = 1:length(All_feature)
%     fprintf(fid,'%.4f %.4f %.4f %.4f %.4f %.4f %.4f %.4f %.4f %.4f %.4f %.4f %.4f\r\n',All_feature(n,1:end));
% end
% fclose(fid);
 save('5.81.txt','All_feature','-ascii');
                                
toc
        
