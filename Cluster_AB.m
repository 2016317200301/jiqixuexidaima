load GSE2034_ma2;
ma2=zscore(ma2);  %%标准化，使用同一个数据集分析，可以不做。
variance=var(ma2')';  %%计算方差
index=(1:length(ma2_geneId))';

index=[index variance];

index=sortrows(index,-2);  %%按照方差从大到小排序

data=ma2(index(1:20,1),:);  %%取出方差最大的20个基因的表达值

patient=anno(:,1);   %%病人ID

gene=ma2_geneName(index(1:20,1),1);  %%用于聚类的20个基因的名字

clustergram(data,'RowLabels',gene,'ColumnLabels',patient);  %%层次聚类
idx=kmeans(data',4);   %%k-means聚类

