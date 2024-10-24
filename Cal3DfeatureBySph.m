clc,clear;
tic;

%%%%%%%%%%%%%%%%%%%�˴�������������Ƶ���ά����(�������򣩣��������������뾶����λΪm)%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%����һ�����ƣ�����õ������������뾶�µ�planarity��ƽ���ԣ���linearity�����ԣ��Լ�sphericity�������ԣ�

Point_data = load ('test.txt'); %��ȡ��������
MD1 = KDTreeSearcher(Point_data); %����KD-treeģ��
[row,line] = size(Point_data);%%��ԭʼ������У�

C_Lambda_all = zeros(row,1);
D_Feature_all = zeros(row,3);
% E_Dim_all = zeros(row,1);

rSize=zeros(10,1);
rSize(1)=0.1;

for i= 2:10
    rSize(i)=0.1+(i^2)*0.03;
end

toc

 for i = 1:row %����ÿ��
     tic
    C_Lambda = zeros(10,1); 
    D_Feature = zeros(10,3);%�����洢�õ㲻ͬ�����LPS
    E_Dim = zeros(10,1);    %�����洢�õ㲻ͬ�������ũ��
  
    for k = 1:10 
        j=Size(k,1);
        Idx = rangesearch(MD1,Point_data(i,:),j);%��Point_data�и�������뾶���ٽ��㣬��������
        NB_point = Point_data(Idx{1}',:);%��Point_data�и��������õ��ٽ��㼯

        [r,l]=size(NB_point);
        S = D_tensor(NB_point,r);  %%����Э������
        [Z,L] = eig(S);          %%����ֵ
        Normal_v(1,:) = (Z(:,1))';   %�����������ھ���Normal_v��

        [d,ind]=sort(diag(L),'descend'); %%������ֵ������ķ�������
        Lambda =d';

% %%%%%%%������ķ�����������ֵȡ����%%%%%%%%%%%%%%
%          e1=sqrt(Lambda(1));
%          e2=sqrt(Lambda(2));
%          e3=sqrt(Lambda(3));
        
                
%%%%%%%%%%%%%%%%%%%% ԭʼ���룬��������ֵ��һ��%%%%%%%%%%%%%%
        lambda_sum = sum(Lambda);
        e1 = Lambda(1)/lambda_sum;
        e2 = Lambda(2)/lambda_sum;
        e3 = Lambda(3)/lambda_sum;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        L_lambda = (e1-e2) / e1;
        P_lambda = (e2-e3) /e1;
        S_lambda = e3  /e1 ; %����ά������L_lambda��P_lambda��S_lambda
        Feature = [L_lambda,P_lambda,S_lambda]; 
        
        E_dim = -L_lambda*log(L_lambda)-P_lambda*log(P_lambda)-S_lambda*log(S_lambda);
        D_Feature(k,:) = Feature;
        E_Dim(k,1) = E_dim;
    end
        [~,OT] = min(E_Dim);   %�ҵ���С����ũ�ض�Ӧ���±�
        D_Feature_all(i,:) = D_Feature(OT,:);
%       E_Dim_all(i,:) = E_Dim(OT,:);
%         
  if mod(i,1000) == 0 
        fprintf('��%d����\n',i);
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
        
