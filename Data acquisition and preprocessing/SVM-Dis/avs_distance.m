%�����о����D,y=Et.����cavv
%������˹�˲�
function [avs] = avs_distance(inputA,inputB,inputC,output2,dt,T,JJ)
%UNTITLED Summary of this function goes here
%Detailed explanation goes here
%inputA ����UD���ٶ�
%inputB ����EW���ٶ�
%inputC ����NS���ٶ�
%output2 P����ʱ
%dt ʱ����
[m n]=size(inputA);%������ٶȵľ���ά��
for j=1:n
    t(j,1)=(j-1)*dt;
end
t4=T;%��ʾʱ�䴰����s
v=cumtrapz(t,inputA);%�ٶȣ����Ƕ�UD���ٶ������
v1=cumtrapz(t,inputB);%�ٶȣ����Ƕ�EW���ٶ������
v2=cumtrapz(t,inputC);%�ٶȣ����Ƕ�NS���ٶ������
%d=cumtrapz(t,v);%λ�ƣ����Ƕ��ٶ������
%----------------4��0.075HZ-3HZ��ͨ�˲�-----------------------------------------%
                                                      %�˲��׶�
fmax=45;      fmin=40;
cutup=fmax/(1/(2*dt));
cutdown=fmin/(1/(2*dt));
wn=[cutdown, cutup]; 
[b,a]=butter(JJ,wn);                                               %4���˲�
%----------------4��0.075HZ-3HZ��ͨ�˲�-----------------------------------------%


%----------------4��0.075HZ��ͨ�˲�-----------------------------------------%
% Fs = 1/dt;
% fl= 0.075; % low cutoff frequency for highpass filter, depending on you
% ftype = 'high';
% f_order = 4;% order
% [b,a] = butter(f_order,fl/(Fs/2),ftype);% f_order order butterworth highpass filter     
%----------------4��0.075HZ��ͨ�˲�-----------------------------------------%
v=filter(b,a,v);%�˲����UD�ٶ�
v1=filter(b,a,v1);%�˲����EW�ٶ�
v2=filter(b,a,v2);%�˲����NS�ٶ�



[~,q]=min(abs(t-output2));
n1=find(t==t(q,1));%n1��P����ʱ�ĺ�һ����
t2=0:dt:t4;
n1=n1-2;

k=t4/dt;
caaA=abs(v(n1:n1+k-1));%caa�ۻ����Լ��ٶ�UD
caaB=abs(v1(n1:n1+k-1));%caa�ۻ����Լ��ٶ�UD
caaC=abs(v2(n1:n1+k-1));%caa�ۻ����Լ��ٶ�UD
hecheng=sqrt(((caaA).^2)+((caaB).^2)+((caaC).^2));%ʱ�䴰��������ϳɼ��ٶ�
avs=sum(abs(hecheng));%avs�ۻ�����λ��UD


end

