%����׿Խ���ڣ���ֵλ�ƣ���ֵ�ٶȣ���ֵ���ٶ�
%������˹�˲�
function [Pa,Pv]= Pa_Pv_lvbo_wave(inputA,inputB,inputC,output2,dt)
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
inputB1=filter(b,a,inputB);%�˲���ļ��ٶ�
inputC1=filter(b,a,inputC);%�˲���ļ��ٶ�

vA=filter(b,a,vA);%�˲�����ٶ�
vB=filter(b,a,vB);%�˲�����ٶ�
vC=filter(b,a,vC);%�˲�����ٶ�




[~,q1]=min(abs(t1-output2));
[~,q2]=min(abs(t2-output2));
[~,q3]=min(abs(t3-output2));
% t1=output2:output2+3;
nn1=find(t1==t1(q1,1));%n1��P����ʱ�ĺ�һ����
nn2=find(t2==t2(q2,1));%n1��P����ʱ�ĺ�һ����
nn3=find(t3==t3(q3,1));%n1��P����ʱ�ĺ�һ����


%%%��ֵ�ٶȵĲ���%%%%%
k=t4/dt;
for i=1:1:k
    vud=abs(vA(nn1:nn1+i));
    vew=abs(vB(nn2:nn2+i));
    vns=abs(vC(nn3:nn3+i));
    Pv(i)= max(sqrt(((vud).^2)+((vew).^2)+((vns).^2)));%ʱ�䴰��������ϳɼ��ٶ�
end
%%%��ֵ�ٶȵĲ���%%%%%

%%%��ֵ���ٶȵĲ���%%%%%
k=t4/dt;
for i=1:1:k
    aud=abs(inputA1(nn1:nn1+i));
    aew=abs(inputB1(nn2:nn2+i));
    ans=abs(inputC1(nn3:nn3+i));
    Pa(i)= max(sqrt(((aud).^2)+((aew).^2)+((ans).^2)));%ʱ�䴰��������ϳɼ��ٶ�
end
%%%��ֵ���ٶȵĲ���%%%%%




end

