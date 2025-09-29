%�����о����A,y=At. ���ݰ�����˹�Ҷ�
%������˹�˲�
function [Ia]= Ia_distance(inputA,inputB,inputC,output2,dt,T,J)
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
v=cumtrapz(t,inputA);%�ٶȣ����ǶԼ��ٶ������
d=cumtrapz(t,v);%λ�ƣ����Ƕ��ٶ������
%----------------4��0.075HZ-3HZ��ͨ�˲�-----------------------------------------%
                                                      %�˲��׶�
fmax=45;      fmin=40;
cutup=fmax/(1/(2*dt));
cutdown=fmin/(1/(2*dt));
wn=[cutdown, cutup]; 
[b,a]=butter(J,wn);                                               %4���˲�
%----------------4��0.075HZ-3HZ��ͨ�˲�-----------------------------------------%


%----------------4��0.075HZ��ͨ�˲�-----------------------------------------%
% Fs = 1/dt;
% fl= 0.075; % low cutoff frequency for highpass filter, depending on you
% ftype = 'high';
% f_order = 4;% order
% [b,a] = butter(f_order,fl/(Fs/2),ftype);% f_order order butterworth highpass filter     
%----------------4��0.075HZ��ͨ�˲�-----------------------------------------%
inputA1=filter(b,a,inputA);%�˲����UD���ٶ�
inputB1=filter(b,a,inputB);%�˲����EW���ٶ�
inputC1=filter(b,a,inputC);%�˲����NS���ٶ�

[~,q]=min(abs(t-output2));

n1=find(t==t(q,1));%n1��P����ʱ�ĺ�һ����
t2=0:dt:t4;
n1=n1-2;

inputA1=inputA1(n1:n1+t4/dt);%ʱ�䴰�ڵ�UD���ٶ�
inputB1=inputB1(n1:n1+t4/dt);%ʱ�䴰�ڵ�EW���ٶ�
inputC1=inputC1(n1:n1+t4/dt);%ʱ�䴰�ڵ�NS���ٶ�

syn_acc= sqrt(((inputA1).^2)+((inputB1).^2)+((inputC1).^2));%������ϳɼ��ٶ�
y=syn_acc;
Ia=(pi/(2*9.8*100))*max((cumtrapz(t2,y.^2)));%������˹�Ҷ�Ia


end

