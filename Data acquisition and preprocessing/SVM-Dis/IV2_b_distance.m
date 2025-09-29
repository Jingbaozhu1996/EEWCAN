%�����о����D,y=Et.����IV2
%������˹�˲�
function [IV2_b]= IV2_b_distance(inputA,inputB,inputC,output2,dt,T,J)
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
%inputA1=filter(b,a,inputA);%�˲����UD���ٶ�
v=filter(b,a,v);%�˲����UD�ٶ�
v1=filter(b,a,v1);%�˲����UD�ٶ�
v2=filter(b,a,v2);%�˲����UD�ٶ�
syn_vel= sqrt(((v).^2)+((v1).^2)+((v2).^2));%������ϳɼ��ٶ�

yy=cumtrapz(t,syn_vel.^2);%�ٶ�ƽ�����֣��������

% y=syn_vel+0.001;
% 
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

%a1=y(n1:n1+t4/dt);%ʱ�䴰�ڵ��ٶ�
MM=yy(n1-1:n1-2+t4/dt);

YY=[];
tt=0.01;
ttt=[];
a=0;

for i=1:t4/dt
    %MM(i)=max(y(n1:n1+i-1));  %%����������ϵĵ�
    YY(i)=log(MM(i)/tt);    
    tt=tt+0.01;
    ttt(i)=a+0.01;
    a=a+0.01;
end

ttt=ttt';
YY=YY';
X=[ones(length(ttt),1),ttt];Y=YY; %y=lnB+At
[b,bint,r,rint,stats]=regress(Y,X); 
IV2_b=exp(b(1));

end

