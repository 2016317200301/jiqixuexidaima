load GSE2034_ma2;
ma2=zscore(ma2);  %%��׼����ʹ��ͬһ�����ݼ����������Բ�����
variance=var(ma2')';  %%���㷽��
index=(1:length(ma2_geneId))';

index=[index variance];

index=sortrows(index,-2);  %%���շ���Ӵ�С����

data=ma2(index(1:20,1),:);  %%ȡ����������20������ı��ֵ

patient=anno(:,1);   %%����ID

gene=ma2_geneName(index(1:20,1),1);  %%���ھ����20�����������

clustergram(data,'RowLabels',gene,'ColumnLabels',patient);  %%��ξ���
idx=kmeans(data',4);   %%k-means����

