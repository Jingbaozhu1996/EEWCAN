%����MMI�Ҷ�PGA\PGV����
%������˹�˲�
function [MMI_ACC,MMI_VEL]= MMI_Jisuan_ACC_VEL(inputA,dt)
%UNTITLED Summary of this function goes here
%Detailed explanation goes here
%inputA ������ٶ�
%output2 P����ʱ
%dt ʱ����
[~, n]=size(inputA);%������ٶȵľ���ά��
for j=1:n
    t(j,1)=(j-1)*dt;
end

v=cumtrapz(t,inputA);%�ٶȣ����ǶԼ��ٶ������

%----------------4��0.075HZ-3HZ��ͨ�˲�-----------------------------------------%
                                                      %�˲��׶�
% fmax=0.1;      fmin=10;
% cutup=fmax/(1/(2*dt));
% cutdown=fmin/(1/(2*dt));
% wn=[cutdown, cutup]; 
% [b,a]=butter(4,wn);                                               %4���˲�
%----------------4��0.075HZ-3HZ��ͨ�˲�-----------------------------------------%


% %----------------4��0.075HZ��ͨ�˲�-----------------------------------------%
Fs = 1/dt;
fl= 0.075; % low cutoff frequency for highpass filter, depending on you
ftype = 'high';
f_order = 4;% order
[b,a] = butter(f_order,fl/(Fs/2),ftype);% f_order order butterworth highpass filter     
% %----------------4��0.075HZ��ͨ�˲�-----------------------------------------%

inputA1=filter(b,a,inputA);%�˲���ļ��ٶ�
MMI_ACC=inputA1;%��ֵ���ٶ�

v=filter(b,a,v);%�˲�����ٶ�
MMI_VEL=v;%��ֵ�ٶ�



end

