%����׿Խ���ڣ���ֵλ�ƣ���ֵ�ٶȣ���ֵ���ٶ�
%������˹�˲�
function [Tvaud]= Tvaud_lvbo_wave(inputA,output2,dt)
%UNTITLED Summary of this function goes here  
%   Detailed explanation goes here
%inputA ������ٶ�
%output2 P����ʱ
%dt ʱ����
[m n1]=size(inputA);%������ٶȵľ���ά��
for j=1:n1
    t1(j,1)=(j-1)*dt;
end



t4=3;%��ʾʱ�䴰����s
vA=cumtrapz(t1,inputA);%�ٶȣ����ǶԼ��ٶ������
dA=cumtrapz(t1,vA);%λ�ƣ����Ƕ��ٶ������


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

vA=filter(b,a,vA);%�˲�����ٶ�

dA=filter(b,a,dA);%�˲����λ��


[~,q1]=min(abs(t1-output2));

% t1=output2:output2+3;
nn1=find(t1==t1(q1,1));%n1��P����ʱ�ĺ�һ����


t2=0:dt:t4;



%%%��ֵ��Tva%%%%%%%%%%%%%%%%%%%%% 
k=t4/dt;
for i=1:1:k
    vud=abs(vA(nn1:nn1+i));
    
    aud=abs(inputA1(nn1:nn1+i));

    Tvaud(i)=2*pi*(vud/aud);%P����ʱ��ʱ�䴰�ڵķ�ֵ��
end
%%%��ֵ��Tva%%%%%%%%%%%%%%%%%%%%%  


end

