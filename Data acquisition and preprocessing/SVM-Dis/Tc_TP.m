%����׿Խ���ڣ���ֵλ�ƣ���ֵ�ٶȣ���ֵ���ٶ�
%������˹�˲�
function [tc,TP]= Tc_TP(inputA,output2,dt,T)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%inputA ������ٶ�
%output2 P����ʱ
%dt ʱ����
[m n]=size(inputA);%������ٶȵľ���ά��
for j=1:n
    t(j,1)=(j-1)*dt;
end
t4=T;%��ʾʱ�䴰����s
v=cumtrapz(t,inputA);%�ٶȣ����ǶԼ��ٶ������
d=cumtrapz(t,v);%λ�ƣ����Ƕ��ٶ������
%----------------4��0.075HZ-3HZ��ͨ�˲�-----------------------------------------%
                                                      %�˲��׶�
% fmax=3;      fmin=0.075;
% cutup=fmax/(1/(2*dt));
% cutdown=fmin/(1/(2*dt));
% wn=[cutdown, cutup]; 
% [b,a]=butter(4,wn);                                               %4���˲�
%----------------4��0.075HZ-3HZ��ͨ�˲�-----------------------------------------%


%----------------4��0.075HZ��ͨ�˲�-----------------------------------------%
Fs = 1/dt;
fl= 0.075; % low cutoff frequency for highpass filter, depending on you
ftype = 'high';
f_order = 4;% order
[b,a] = butter(f_order,fl/(Fs/2),ftype);% f_order order butterworth highpass filter     
%----------------4��0.075HZ��ͨ�˲�-----------------------------------------%
inputA1=filter(b,a,inputA);%�˲���ļ��ٶ�
d=filter(b,a,d);%�˲����λ��
v=filter(b,a,v);%�˲�����ٶ�
[~,q]=min(abs(t-output2));
% t1=output2:output2+3;
n1=find(t==t(q,1));%n1��P����ʱ�ĺ�һ����
% v1=v(n1:);
% d1=d(n1:);
% a1=inputA1(n1:n1+t4/0.01);
t2=0:dt:t4;
% v1=cumtrapz(t2,a1);
% d1=cumtrapz(t2,v1);
if n1+t4/dt<=n
    v1=v(n1:n1+t4/dt);%ʱ�䴰�ڵ��ٶ�
    d1=d(n1:n1+t4/dt);%ʱ�䴰�ڵ�λ��
    a1=inputA1(n1:n1+t4/dt);%ʱ�䴰�ڵļ��ٶ�
else
    v1=v(n1:n);%ʱ�䴰�ڵ��ٶ�
    d1=d(n1:n);%ʱ�䴰�ڵ�λ��
    a1=inputA1(n1:n);%ʱ�䴰�ڵļ��ٶ�
    t2=0:dt:(n-n1)*dt;
end

Pd=max(abs(d1));%P����ʱ��ʱ�䴰�ڵķ�ֵλ��
% r1=quad(inline('(v.^2)'),0,4);
% r2=quad(inline('(d.^2)'),0,4);
r1=trapz(t2,v1.^2);
r2=trapz(t2,d1.^2);
r=r1/r2;
tc=2*pi/sqrt(r);%��������
TP=tc*Pd;%�������

    
end

