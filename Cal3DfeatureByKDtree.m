clc,clear;
tic;

%%%%%%%%%%%%%%%%%%%�˴�������������Ƶ���ά�������������������뾶���ٽ�������%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%����һ�����ƣ�����õ������������뾶�µ�planarity��ƽ���ԣ���linearity�����ԣ��Լ�sphericity�������ԣ�


Point_data =load('merge.txt'); %��ȡ��������

MD1 = KDTreeSearcher(Point_data); %����KD-treeģ��

[row,line] = size(Point_data);%%��ԭʼ������У�

%%�ҵ��ٽ���150���㣬����IDX�洢������Щ����Point_data�е��кţ�������Ϊ�����ڵĵ�ţ�,
%%D�洢���Ǿ��룬������Ϊ���������ľ��룬�������У�����Խ��ԽԶ
%%����knnsearch�ҵ�����Ŀ��������150���ٽ���
[IDX,D] = knnsearch(MD1,Point_data,'K',150); 


C_Lambda_all = zeros(row,1);
D_Feature_all = zeros(row,3);
E_Dim_all = zeros(row,1);


 for i = 1:row %����ÿ��
    
    D_Feature = zeros(10,3);%�����洢�õ㲻ͬ�����ά������
    E_Dim = zeros(10,1);    %�����洢�õ㲻ͬ�������ũ��
    E_Lambda = zeros(10,1);
         
    for j = 15:15:150 %%������Ϊ10��������ũ�ؼ�����ٸ����������ѣ�
        ord=j/15;
        NB_point = Point_data(IDX(i,1:j),:);%��Point_data�и���IDX�ҵ������i���������K���ٽ�������겢�þ���NB_point�洢
%IDX����Ϊ��Ϊ��ţ���Ϊ�õ��100���ٽ����


        S = D_tensor(NB_point,j);  %%����Э���������ú���D_tensor���i�����K���ٽ������3D�ṹ����S
        [Z,L] = eig(S);            %%����ֵ
        Normal_v(1,:) = (Z(:,1))';   %�����������ھ���Normal_v��

        [d,ind]=sort(diag(L),'descend'); %%������ֵ������ķ�������
%       Lambda = [L(3,3),L(2,2),L(1,1)]; %�����i�������Ӧ��K���ٽ����3D�ṹ����S������ֵ���þ���Lambda�洢
        Lambda =d';

  

%%%%%%%%%%%%%%%%%%%% ����ֵ��һ��%%%%%%%%%%%%%%
        lambda_sum = sum(Lambda);
        e1 = Lambda(1)/lambda_sum;
        e2 = Lambda(2)/lambda_sum;
        e3 = Lambda(3)/lambda_sum;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        L_lambda = (e1-e2) / e1;
        P_lambda = (e2-e3) /e1;
        S_lambda = e3  /e1 ; %����ά������L��P��S
        
        Feature = [L_lambda,P_lambda,S_lambda]; 
        E_dim = -L_lambda*log(L_lambda)-P_lambda*log(P_lambda)-S_lambda*log(S_lambda);
        
        D_Feature(ord,:) = Feature;
        E_Dim(ord,1) = E_dim;
    end
        [~,OT] = min(E_Dim);   %�ҵ���i����K��10ȡ��100���õ�����С����ũ�ض�Ӧ���±�
        D_Feature_all(i,:) = D_Feature(OT,:);
        E_Dim_all(i,:) = E_Dim(OT,:);
        
   if mod(i,1000) == 0 
         fprintf('��%d����\n',i);

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

