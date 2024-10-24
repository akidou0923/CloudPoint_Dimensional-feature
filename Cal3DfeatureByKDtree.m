clc,clear;
tic;

%%%%%%%%%%%%%%%%%%%此代码用来计算点云的三维特征，并给出最佳邻域半径（临近点数）%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%输入一个点云，输出该点云在最佳邻域半径下的planarity（平面性）、linearity（线性）以及sphericity（立体性）


Point_data =load('merge.txt'); %读取点云数据

MD1 = KDTreeSearcher(Point_data); %建立KD-tree模型

[row,line] = size(Point_data);%%（原始点的行列）

%%找的临近的150个点，其中IDX存储的是这些点在Point_data中的行号（行向量为邻域内的点号）,
%%D存储的是距离，行向量为点与邻域点的距离，升序排列，距离越来越远
%%利用knnsearch找到距离目标点最近的150个临近点
[IDX,D] = knnsearch(MD1,Point_data,'K',150); 


C_Lambda_all = zeros(row,1);
D_Feature_all = zeros(row,3);
E_Dim_all = zeros(row,1);


 for i = 1:row %遍历每点
    
    D_Feature = zeros(10,3);%用来存储该点不同邻域的维度特征
    E_Dim = zeros(10,1);    %用来存储该点不同邻域的香农熵
    E_Lambda = zeros(10,1);
         
    for j = 15:15:150 %%（步长为10，根据香农熵计算多少个点的邻域最佳）
        ord=j/15;
        NB_point = Point_data(IDX(i,1:j),:);%从Point_data中根据IDX找到距离第i个点最近的K个临近点的坐标并用矩阵NB_point存储
%IDX排列为行为点号，列为该点的100个临近点号


        S = D_tensor(NB_point,j);  %%计算协方差阵，利用函数D_tensor与第i个点的K个临近点计算3D结构张量S
        [Z,L] = eig(S);            %%特征值
        Normal_v(1,:) = (Z(:,1))';   %将法向量存在矩阵Normal_v中

        [d,ind]=sort(diag(L),'descend'); %%将特征值按降序的方法排列
%       Lambda = [L(3,3),L(2,2),L(1,1)]; %计算第i个点与对应的K个临近点的3D结构张量S的特征值并用矩阵Lambda存储
        Lambda =d';

  

%%%%%%%%%%%%%%%%%%%% 特征值归一化%%%%%%%%%%%%%%
        lambda_sum = sum(Lambda);
        e1 = Lambda(1)/lambda_sum;
        e2 = Lambda(2)/lambda_sum;
        e3 = Lambda(3)/lambda_sum;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        L_lambda = (e1-e2) / e1;
        P_lambda = (e2-e3) /e1;
        S_lambda = e3  /e1 ; %计算维度特征L、P、S
        
        Feature = [L_lambda,P_lambda,S_lambda]; 
        E_dim = -L_lambda*log(L_lambda)-P_lambda*log(P_lambda)-S_lambda*log(S_lambda);
        
        D_Feature(ord,:) = Feature;
        E_Dim(ord,1) = E_dim;
    end
        [~,OT] = min(E_Dim);   %找到第i个点K从10取到100所得到的最小的香农熵对应的下标
        D_Feature_all(i,:) = D_Feature(OT,:);
        E_Dim_all(i,:) = E_Dim(OT,:);
        
   if mod(i,1000) == 0 
         fprintf('第%d个点\n',i);

  end
 
    
 end
 
All_feature = [Point_data(:,:),D_Feature_all];
% fid = fopen('feature_3D.txt','wt');
% for n = 1:length(All_feature)
%     fprintf(fid,'%.4f %.4f %.4f %.4f %.4f %.4f %.4f %.4f %.4f %.4f %.4f %.4f %.4f\r\n',All_feature(n,1:end));
% end
% fclose(fid);
save('result.txt','All_feature','-ascii');
                                
toc

