load GSE2034_ma2;
data=ma2(:,co>=0);  
label=co(co>=0,1);
N=length(data(1,:));

attr_num=length(data(:,1));
index=(1:attr_num)';


Accuracy_svm=0;  %初始化变量，存储分类器的准确率
auc_svm=0;       %初始化变量，存储分类器的auc
sensitivity_svm=0; %初始化变量，存储分类器的灵敏度
specificity_svm=0; %初始化变量，存储分类器的特异度


for i=1:10
    Indices = crossvalind('Kfold', N, 5);   %%五倍交叉验证，将样本分成五份
    for j=1:5
        train=data(:,Indices~=j);     %%训练集数据
        train_label=label(Indices~=j,1);  %%训练集标签
        test=data(:,Indices==j);     %%测试集数据
        test_label=label(Indices==j,1);  %%测试集标签
        
        [p,t]=mattest(train(:,train_label==1),train(:,train_label==0));   %%在训练集中用t-test选择差异表达基因
        feature_loc=[index p];
        feature_loc=sortrows(feature_loc,2);
        feature=feature_loc(1:100,1); %假设选择差异最大的100个属性作为特征
        train_data=(train(feature,:))';
        test_data=(test(feature,:))';
        %%SVM
        nb_svm=fitcsvm(train_data,train_label);    %%SVM的训练模型
        [res_svm, score] = predict(nb_svm,test_data);          %%SVM模型去预测病人的类别，res_SVM是标签值，score是一个连续的score，也就是每个样本属于类别1的得分
        if var(score)~=0
            [~,~,~,auc]=perfcurve(test_label,res_svm,'1'); %计算分类器的AUC
        else
            auc=0.5;
        end
        auc_svm=auc_svm+auc;
        cpre=classperf(test_label,res_svm);  %计算分类器的准确率等指标
        Accuracy_svm=Accuracy_svm+cpre.CorrectRate;  %%准确率
        sensitivity_svm=sensitivity_svm+cpre.Sensitivity;  %%灵敏度
        specificity_svm=specificity_svm+cpre.Specificity;   %%特异度
    end
end

result=[auc_svm Accuracy_svm sensitivity_svm specificity_svm]./50;

title={'AUC','Accuracy','Seneitivity','Specificity'};
res=[title;num2cell(result)];

xlswrite('Classification_result.xlsx',res);


