
function [a1,d1,v1]= Dataa(inputA,output2,dt)

[m n]=size(inputA);%������ٶȵľ���ά��
for j=1:n
    t(j,1)=(j-1)*dt;
end
t4=3;%��ʾʱ�䴰����s
v=cumtrapz(t,inputA);%�ٶȣ����ǶԼ��ٶ������
d=cumtrapz(t,v);%λ�ƣ����Ƕ��ٶ������
[~,q]=min(abs(t-output2));
% t1=output2:output2+3;
n1=find(t==t(q,1));%n1��P����ʱ�ĺ�һ����
v1=abs(v(n1+1:n1+t4/dt));%p����ʱ��ʱ�䴰�ڵ��ٶȷ�ֵ
d1=abs(d(n1+1:n1+t4/dt));%ʱ�䴰�ڵ�λ��
a1=abs(inputA(n1+1:n1+t4/dt));%P����3s�ļ��ٶ�
if n1-2-t4/dt<=0
    v2=max(abs(v(1:n1-2)));%P��ǰ3s���ٶȷ�ֵ
    a2=inputA(1:n1-2);
else
    v2=max(abs(v(n1-2-t4/dt:n1-2)));%P��ǰ3s���ٶȷ�ֵ
    a2=inputA(n1-2-t4/dt:n1-2);
end
SNR=v1/v2;%�����
end

