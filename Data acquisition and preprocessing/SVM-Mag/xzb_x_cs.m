%����׿Խ���ڣ���ֵλ�ƣ���ֵ�ٶȣ���ֵ���ٶ�
%������˹�˲�
function [Pv,Pa,IV2,Tva,cav,caa,DE]= xzb_x_cs(inputA,output2,dt)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%inputA ������ٶ�
%output2 P����ʱ
%dt ʱ����
[m n]=size(inputA);%������ٶȵľ���ά��
for j=1:n
    t(j,1)=(j-1)*dt;
end
t4=3;%��ʾʱ�䴰����s
v=cumtrapz(t,inputA);%�ٶȣ����ǶԼ��ٶ������
d=cumtrapz(t,v);%λ�ƣ����Ƕ��ٶ������
%----------------4��0.075HZ-3HZ��ͨ�˲�-----------------------------------------%
                                                      %�˲��׶�
fmax=3;      fmin=0.075;
cutup=fmax/(1/(2*dt));
cutdown=fmin/(1/(2*dt));
wn=[cutdown, cutup]; 
[b,a]=butter(4,wn);                                               %4���˲�
%----------------4��0.075HZ-3HZ��ͨ�˲�-----------------------------------------%


% %----------------4��0.075HZ��ͨ�˲�-----------------------------------------%
% Fs = 1/dt;
% fl= 0.075; % low cutoff frequency for highpass filter, depending on you
% ftype = 'high';
% f_order = 4;% order
% [b,a] = butter(f_order,fl/(Fs/2),ftype);% f_order order butterworth highpass filter     
% %----------------4��0.075HZ��ͨ�˲�-----------------------------------------%

inputA1=filter(b,a,inputA);%�˲���ļ��ٶ�

v=filter(b,a,v);%�˲�����ٶ�
d=filter(b,a,d);%�˲����λ��
[~,q]=min(abs(t-output2));
% t1=output2:output2+3;
n1=find(t==t(q,1));%n1��P����ʱ�ĺ�һ����
% v1=v(n1:);
% d1=d(n1:);
% a1=inputA1(n1:n1+t4/0.01);
t2=0:dt:t4;

if n1+t4/dt<=n
    v1=v(n1:n1+t4/dt);%ʱ�䴰�ڵ��ٶ�
    d1=d(n1:n1+t4/dt);%ʱ�䴰�ڵ�λ��
    a1=inputA1(n1:n1+t4/dt);%ʱ�䴰�ڵļ��ٶ�
    t3=output2:dt:(output2+t4);
    IV2=trapz(t3,v1.^2);%�ٶ�ƽ������
else
    v1=v(n1:n);%ʱ�䴰�ڵ��ٶ�
    d1=d(n1:n);%ʱ�䴰�ڵ�λ��
    a1=inputA1(n1:n);%ʱ�䴰�ڵļ��ٶ� 
    t2=0:dt:(n-n1)*dt;
    t3=output2:dt:(output2+(n-n1)*dt);
    IV2=trapz(t3,v1.^2);%�ٶ�ƽ������
end

Pv=max(abs(v1));%P����ʱ��ʱ�䴰�ڵķ�ֵ�ٶ�
Pa=max(abs(a1));%P����ʱ��ʱ�䴰�ڵķ�ֵ���ٶ�
Tva=2*pi*(Pv/Pa);%P����ʱ��ʱ�䴰�ڵķ�ֵ��
%--------------�����ٶ�ƽ������IV2--------------------------%

%--------------�����ٶ�ƽ������IV2--------------------------%
%������һ�������ϵ�
cav=sum(abs(v1));%CAV�ۻ������ٶ�UD
caa=sum(abs(a1));%CAA�ۻ����Լ��ٶ�UD
DI=log(abs(a1.*v1));
DE=max(DI);%�ۻ������仯��


end

