%计算卓越周期，峰值位移，峰值速度，峰值加速度
%巴特沃斯滤波
function [v1,a1]= raw_lvbo_acc_vel(inputA,inputB,inputC,output2,dt)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%inputA 输入加速度
%output2 P波到时
%dt 时间间隔
[m n1]=size(inputA);%输入加速度的矩阵维度
for j=1:n1
    t1(j,1)=(j-1)*dt;
end

[m n2]=size(inputB);%输入加速度的矩阵维度
for j=1:n2
    t2(j,1)=(j-1)*dt;
end

[m n3]=size(inputC);%输入加速度的矩阵维度
for j=1:n3
    t3(j,1)=(j-1)*dt;
end
t4=3;%表示时间窗长度s

vA=cumtrapz(t1,inputA);%速度，这是对加速度求积分
dA=cumtrapz(t1,vA);%位移，这是对速度求积分

vB=cumtrapz(t2,inputB);%速度，这是对加速度求积分
dB=cumtrapz(t2,vB);%位移，这是对速度求积分

vC=cumtrapz(t3,inputC);%速度，这是对加速度求积分
dC=cumtrapz(t3,vC);%位移，这是对速度求积分
%----------------4阶0.075HZ-3HZ带通滤波-----------------------------------------%
                                                      %滤波阶段
fmax=3;      fmin=0.075;
cutup=fmax/(1/(2*dt));
cutdown=fmin/(1/(2*dt));
wn=[cutdown, cutup]; 
[b,a]=butter(4,wn);                                               %4阶滤波
%----------------4阶0.075HZ-3HZ带通滤波-----------------------------------------%


%----------------4阶0.075HZ高通滤波-----------------------------------------%
% Fs = 1/dt;
% fl= 0.075; % low cutoff frequency for highpass filter, depending on you
% ftype = 'high';  %low
% f_order = 4;% order
% [b,a] = butter(f_order,fl/(Fs/2),ftype);% f_order order butterworth highpass filter     
%----------------4阶0.075HZ高通滤波-----------------------------------------%
inputA1=filter(b,a,inputA);%滤波后的加速度
inputB1=filter(b,a,inputB);
inputC1=filter(b,a,inputC);

vA=filter(b,a,vA);%滤波后的速度
vB=filter(b,a,vB);%滤波后的速度
vC=filter(b,a,vC);%滤波后的速度

dA=filter(b,a,dA);%滤波后的位移
dB=filter(b,a,dB);%滤波后的位移
dC=filter(b,a,dC);%滤波后的位移

[~,q1]=min(abs(t1-output2));
[~,q2]=min(abs(t2-output2));
[~,q3]=min(abs(t3-output2));
% t1=output2:output2+3;
nn1=find(t1==t1(q1,1));%n1是P波到时的后一个点
nn2=find(t2==t2(q2,1));%n1是P波到时的后一个点
nn3=find(t3==t3(q3,1));%n1是P波到时的后一个点

t2=0:dt:t4;
vud=abs(vA(nn1+1:nn1+t4/dt));%时间窗内的位移
vew=abs(vB(nn2+1:nn2+t4/dt));%时间窗内的位移
vns=abs(vC(nn3+1:nn3+t4/dt));%时间窗内的位移
v1= sqrt(((vud).^2)+((vew).^2)+((vns).^2));%时间窗内三分向合成位移

aud=abs(inputA1(nn1+1:nn1+t4/dt));%时间窗内的位移
aew=abs(inputB1(nn2+1:nn2+t4/dt));%时间窗内的位移
ans=abs(inputC1(nn3+1:nn3+t4/dt));%时间窗内的位移
a1= sqrt(((aud).^2)+((aew).^2)+((ans).^2));%时间窗内三分向合成位移



end

