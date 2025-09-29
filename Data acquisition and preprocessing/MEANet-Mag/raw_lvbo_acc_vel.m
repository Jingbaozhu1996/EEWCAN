%����׿Խ���ڣ���ֵλ�ƣ���ֵ�ٶȣ���ֵ���ٶ�
%������˹�˲�
function [v1,a1]= raw_lvbo_acc_vel(inputA,inputB,inputC,output2,dt)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%inputA ������ٶ�
%output2 P����ʱ
%dt ʱ����
[m n1]=size(inputA);%������ٶȵľ���ά��
for j=1:n1
    t1(j,1)=(j-1)*dt;
end

[m n2]=size(inputB);%������ٶȵľ���ά��
for j=1:n2
    t2(j,1)=(j-1)*dt;
end

[m n3]=size(inputC);%������ٶȵľ���ά��
for j=1:n3
    t3(j,1)=(j-1)*dt;
end
t4=3;%��ʾʱ�䴰����s

vA=cumtrapz(t1,inputA);%�ٶȣ����ǶԼ��ٶ������
dA=cumtrapz(t1,vA);%λ�ƣ����Ƕ��ٶ������

vB=cumtrapz(t2,inputB);%�ٶȣ����ǶԼ��ٶ������
dB=cumtrapz(t2,vB);%λ�ƣ����Ƕ��ٶ������

vC=cumtrapz(t3,inputC);%�ٶȣ����ǶԼ��ٶ������
dC=cumtrapz(t3,vC);%λ�ƣ����Ƕ��ٶ������
%----------------4��0.075HZ-3HZ��ͨ�˲�-----------------------------------------%
                                                      %�˲��׶�
fmax=3;      fmin=0.075;
cutup=fmax/(1/(2*dt));
cutdown=fmin/(1/(2*dt));
wn=[cutdown, cutup]; 
[b,a]=butter(4,wn);                                               %4���˲�
%----------------4��0.075HZ-3HZ��ͨ�˲�-----------------------------------------%


%----------------4��0.075HZ��ͨ�˲�-----------------------------------------%
% Fs = 1/dt;
% fl= 0.075; % low cutoff frequency for highpass filter, depending on you
% ftype = 'high';  %low
% f_order = 4;% order
% [b,a] = butter(f_order,fl/(Fs/2),ftype);% f_order order butterworth highpass filter     
%----------------4��0.075HZ��ͨ�˲�-----------------------------------------%
inputA1=filter(b,a,inputA);%�˲���ļ��ٶ�
inputB1=filter(b,a,inputB);
inputC1=filter(b,a,inputC);

vA=filter(b,a,vA);%�˲�����ٶ�
vB=filter(b,a,vB);%�˲�����ٶ�
vC=filter(b,a,vC);%�˲�����ٶ�

dA=filter(b,a,dA);%�˲����λ��
dB=filter(b,a,dB);%�˲����λ��
dC=filter(b,a,dC);%�˲����λ��

[~,q1]=min(abs(t1-output2));
[~,q2]=min(abs(t2-output2));
[~,q3]=min(abs(t3-output2));
% t1=output2:output2+3;
nn1=find(t1==t1(q1,1));%n1��P����ʱ�ĺ�һ����
nn2=find(t2==t2(q2,1));%n1��P����ʱ�ĺ�һ����
nn3=find(t3==t3(q3,1));%n1��P����ʱ�ĺ�һ����

t2=0:dt:t4;
vud=abs(vA(nn1+1:nn1+t4/dt));%ʱ�䴰�ڵ�λ��
vew=abs(vB(nn2+1:nn2+t4/dt));%ʱ�䴰�ڵ�λ��
vns=abs(vC(nn3+1:nn3+t4/dt));%ʱ�䴰�ڵ�λ��
v1= sqrt(((vud).^2)+((vew).^2)+((vns).^2));%ʱ�䴰��������ϳ�λ��

aud=abs(inputA1(nn1+1:nn1+t4/dt));%ʱ�䴰�ڵ�λ��
aew=abs(inputB1(nn2+1:nn2+t4/dt));%ʱ�䴰�ڵ�λ��
ans=abs(inputC1(nn3+1:nn3+t4/dt));%ʱ�䴰�ڵ�λ��
a1= sqrt(((aud).^2)+((aew).^2)+((ans).^2));%ʱ�䴰��������ϳ�λ��



end

