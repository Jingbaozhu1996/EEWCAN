%�����о����D,y=Et.����cavd
%������˹�˲�
function [ads_b]= ads_b_distance(inputA,inputB,inputC,output2,dt,T,JJ)
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

d=cumtrapz(t,v);%�ٶȣ����Ƕ�UD���ٶ������
d1=cumtrapz(t,v1);%�ٶȣ����Ƕ�EW���ٶ������
d2=cumtrapz(t,v2);%�ٶȣ����Ƕ�NS���ٶ������
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
d=filter(b,a,d);%�˲����UD�ٶ�
d1=filter(b,a,d1);%�˲����EW�ٶ�
d2=filter(b,a,d2);%�˲����NS�ٶ�



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

%%%cav,caa%%%%%%%%%%%%%%%%%%%%%
k=t4/dt;
caaA=abs(d(n1:n1+k-1));%caa�ۻ����Լ��ٶ�UD
caaB=abs(d1(n1:n1+k-1));%caa�ۻ����Լ��ٶ�UD
caaC=abs(d2(n1:n1+k-1));%caa�ۻ����Լ��ٶ�UD
hecheng=sqrt(((caaA).^2)+((caaB).^2)+((caaC).^2));%ʱ�䴰��������ϳɼ��ٶ�
yy(1)=0;
for i=1:1:k
    yy(i+1)=hecheng(i)+yy(i);
end
MM=yy(2:k+1);

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
ads_b=exp(b(1));


end

