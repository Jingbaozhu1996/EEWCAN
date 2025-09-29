%求震中距参数D,y=Et.根据cavv
%巴特沃斯滤波
function [avs] = avs_distance(inputA,inputB,inputC,output2,dt,T,JJ)
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
v=filter(b,a,v);%滤波后的UD速度
v1=filter(b,a,v1);%滤波后的EW速度
v2=filter(b,a,v2);%滤波后的NS速度



[~,q]=min(abs(t-output2));
n1=find(t==t(q,1));%n1是P波到时的后一个点
t2=0:dt:t4;
n1=n1-2;

k=t4/dt;
caaA=abs(v(n1:n1+k-1));%caa累积绝对加速度UD
caaB=abs(v1(n1:n1+k-1));%caa累积绝对加速度UD
caaC=abs(v2(n1:n1+k-1));%caa累积绝对加速度UD
hecheng=sqrt(((caaA).^2)+((caaB).^2)+((caaC).^2));%时间窗内三分向合成加速度
avs=sum(abs(hecheng));%avs累积绝对位移UD


end

