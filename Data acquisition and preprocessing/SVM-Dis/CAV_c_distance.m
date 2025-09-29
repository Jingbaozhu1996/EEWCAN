%�����о����D,y=Dt.����CAV
%������˹�˲�
function [CAV_c]= CAV_c_distance(inputA,inputB,inputC,output2,dt,T,J)
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
fmax=40;      fmin=35;
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

syn_acc= sqrt(((inputA1).^2)+((inputB1).^2)+((inputC1).^2));%������ϳɼ��ٶ�

%CAA= sqrt(((a_ud).^2)+((a_ew).^2)+((a_ns).^2));%ʱ�䴰��������ϳɼ��ٶ�
yy=cumtrapz(t,abs(syn_acc));%�ۻ������ٶȣ��������
%y=syn_acc+0.001;

% yy=[];
% yy(1)=y(1);
% for i=2:n
%     yy(i)=yy(i-1)+y(i);
% end

[~,q]=min(abs(t-output2));
% t1=output2:output2+3;
n1=find(t==t(q,1));%n1��P����ʱ�ĺ�һ����
% v1=v(n1:);
% d1=d(n1:);
% a1=inputA1(n1:n1+t4/0.01);
t2=0:dt:t4;
% v1=cumtrapz(t2,a1);
% d1=cumtrapz(t2,v1);
n1=n1-2;
%v1=v(n1:n1+t4/dt);%ʱ�䴰�ڵ��ٶ�
%d1=d(n1:n1+t4/dt);%ʱ�䴰�ڵ�λ��
%a1=y(n1:n1+t4/dt);%ʱ�䴰�ڵļ��ٶ�
MM=yy(n1-1:n1-2+t4/dt);
% for i=1:t4/dt
%     MM(i)=yy(n1+1:n1+i);
% end
ttt=[];
a=0;
for k=1:t4/dt
    ttt(k)=a+0.01;
    a=a+0.01;   
end
ttt=ttt';
MM=MM';
X=[ones(length(ttt),1),ttt];Y=[MM];
[b,bint,r,rint,stats]=regress(Y,X); 
CAV_c=b(2);
end

