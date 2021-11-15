load GSE2034_ma2;
data=ma2(:,co>=0);  
label=co(co>=0,1);
N=length(data(1,:));

attr_num=length(data(:,1));
index=(1:attr_num)';


Accuracy_svm=0;  %��ʼ���������洢��������׼ȷ��
auc_svm=0;       %��ʼ���������洢��������auc
sensitivity_svm=0; %��ʼ���������洢��������������
specificity_svm=0; %��ʼ���������洢�������������


for i=1:10
    Indices = crossvalind('Kfold', N, 5);   %%�屶������֤���������ֳ����
    for j=1:5
        train=data(:,Indices~=j);     %%ѵ��������
        train_label=label(Indices~=j,1);  %%ѵ������ǩ
        test=data(:,Indices==j);     %%���Լ�����
        test_label=label(Indices==j,1);  %%���Լ���ǩ
        
        [p,t]=mattest(train(:,train_label==1),train(:,train_label==0));   %%��ѵ��������t-testѡ����������
        feature_loc=[index p];
        feature_loc=sortrows(feature_loc,2);
        feature=feature_loc(1:100,1); %����ѡ���������100��������Ϊ����
        train_data=(train(feature,:))';
        test_data=(test(feature,:))';
        %%SVM
        nb_svm=fitcsvm(train_data,train_label);    %%SVM��ѵ��ģ��
        [res_svm, score] = predict(nb_svm,test_data);          %%SVMģ��ȥԤ�ⲡ�˵����res_SVM�Ǳ�ǩֵ��score��һ��������score��Ҳ����ÿ�������������1�ĵ÷�
        if var(score)~=0
            [~,~,~,auc]=perfcurve(test_label,res_svm,'1'); %�����������AUC
        else
            auc=0.5;
        end
        auc_svm=auc_svm+auc;
        cpre=classperf(test_label,res_svm);  %�����������׼ȷ�ʵ�ָ��
        Accuracy_svm=Accuracy_svm+cpre.CorrectRate;  %%׼ȷ��
        sensitivity_svm=sensitivity_svm+cpre.Sensitivity;  %%������
        specificity_svm=specificity_svm+cpre.Specificity;   %%�����
    end
end

result=[auc_svm Accuracy_svm sensitivity_svm specificity_svm]./50;

title={'AUC','Accuracy','Seneitivity','Specificity'};
res=[title;num2cell(result)];

xlswrite('Classification_result.xlsx',res);


