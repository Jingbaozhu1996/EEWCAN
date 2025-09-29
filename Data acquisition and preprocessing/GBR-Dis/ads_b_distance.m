%求震中距参数D,y=Et.根据cavd
%巴特沃斯滤波
function [ads_b]= ads_b_distance(inputA,inputB,inputC,output2,dt,T,JJ)
%UNTITLED Summary of this function goes here
%Detailed explanation goes here
%inputA 输入UD加速度
%inputB 输入EW加速度
%inputC 输入NS加速度
%output2 P波到时
%dt 时间间隔
[m n]=size(inputA);%输入加速度的矩阵维度
for j=1:n
    t(j,1)=(j-1)*dt;
end
t4=T;%表示时间窗长度s
v=cumtrapz(t,inputA);%速度，这是对UD加速度求积分
v1=cumtrapz(t,inputB);%速度，这是对EW加速度求积分
v2=cumtrapz(t,inputC);%速度，这是对NS加速度求积分

d=cumtrapz(t,v);%速度，这是对UD加速度求积分
d1=cumtrapz(t,v1);%速度，这是对EW加速度求积分
d2=cumtrapz(t,v2);%速度，这是对NS加速度求积分
%d=cumtrapz(t,v);%位移，这是对速度求积分
%----------------4阶0.075HZ-3HZ带通滤波-----------------------------------------%
                                                      %滤波阶段
fmax=45;      fmin=40;
cutup=fmax/(1/(2*dt));
cutdown=fmin/(1/(2*dt));
wn=[cutdown, cutup]; 
[b,a]=butter(JJ,wn);                                               %4阶滤波
%----------------4阶0.075HZ-3HZ带通滤波-----------------------------------------%


%----------------4阶0.075HZ高通滤波-----------------------------------------%
% Fs = 1/dt;
% fl= 0.075; % low cutoff frequency for highpass filter, depending on you
% ftype = 'high';
% f_order = 4;% order
% [b,a] = butter(f_order,fl/(Fs/2),ftype);% f_order order butterworth highpass filter     
%----------------4阶0.075HZ高通滤波-----------------------------------------%
d=filter(b,a,d);%滤波后的UD速度
d1=filter(b,a,d1);%滤波后的EW速度
d2=filter(b,a,d2);%滤波后的NS速度



[~,q]=min(abs(t-output2));
% t1=output2:output2+3;
n1=find(t==t(q,1));%n1是P波到时的后一个点
% v1=v(n1:);
% d1=d(n1:);
% a1=inputA1(n1:n1+t4/0.01);
t2=0:dt:t4;
% v1=cumtrapz(t2,a1);
% d1=cumtrapz(t2,v1);
n1=n1-2;

%%%cav,caa%%%%%%%%%%%%%%%%%%%%%
k=t4/dt;
caaA=abs(d(n1:n1+k-1));%caa累积绝对加速度UD
caaB=abs(d1(n1:n1+k-1));%caa累积绝对加速度UD
caaC=abs(d2(n1:n1+k-1));%caa累积绝对加速度UD
hecheng=sqrt(((caaA).^2)+((caaB).^2)+((caaC).^2));%时间窗内三分向合成加速度
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
    %MM(i)=max(y(n1:n1+i-1));  %%计算包络线上的点
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

