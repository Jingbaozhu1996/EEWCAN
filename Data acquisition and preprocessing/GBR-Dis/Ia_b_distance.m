%求震中距参数A,y=Bt*exp(-At). 根据阿里亚斯烈度
%巴特沃斯滤波
function [Ia_b]= Ia_b_distance(inputA,inputB,inputC,output2,dt,T,J)
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
v=cumtrapz(t,inputA);%速度，这是对加速度求积分
d=cumtrapz(t,v);%位移，这是对速度求积分
%----------------4阶0.075HZ-3HZ带通滤波-----------------------------------------%
                                                      %滤波阶段
fmax=45;      fmin=40;
cutup=fmax/(1/(2*dt));
cutdown=fmin/(1/(2*dt));
wn=[cutdown, cutup]; 
[b,a]=butter(J,wn);                                               %4阶滤波
%----------------4阶0.075HZ-3HZ带通滤波-----------------------------------------%


%----------------4阶0.075HZ高通滤波-----------------------------------------%
% Fs = 1/dt;
% fl= 0.075; % low cutoff frequency for highpass filter, depending on you
% ftype = 'high';
% f_order = 4;% order
% [b,a] = butter(f_order,fl/(Fs/2),ftype);% f_order order butterworth highpass filter     
%----------------4阶0.075HZ高通滤波-----------------------------------------%
inputA1=filter(b,a,inputA);%滤波后的UD加速度
inputB1=filter(b,a,inputB);%滤波后的EW加速度
inputC1=filter(b,a,inputC);%滤波后的NS加速度

syn_acc= sqrt(((inputA1).^2)+((inputB1).^2)+((inputC1).^2));%三分向合成加速度

y=syn_acc;
yy=(pi/(2*9.8*100))*(cumtrapz(t,y.^2));%阿里亚斯烈度Ia

[~,q]=min(abs(t-output2));

n1=find(t==t(q,1));%n1是P波到时的后一个点

t2=0:dt:t4;

n1=n1-2;
v1=v(n1:n1+t4/dt);%时间窗内的速度
d1=d(n1:n1+t4/dt);%时间窗内的位移
a1=y(n1:n1+t4/dt);%时间窗内的加速度

MM=yy(n1:n1-1+t4/dt);

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
Ia_b=exp(b(1));

end

