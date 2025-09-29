%计算卓越周期，峰值位移，峰值速度，峰值加速度
%巴特沃斯滤波
function [omg_ud]= omgud_lvbo_wave(inputA,inputB,inputC,output2,dt)
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
inputB1=filter(b,a,inputB);%滤波后的加速度
inputC1=filter(b,a,inputC);%滤波后的加速度


vA=filter(b,a,vA);%滤波后的速度
vB=filter(b,a,vB);%滤波后的速度
vC=filter(b,a,vC);%滤波后的速度


[~,q1]=min(abs(t1-output2));
[~,q2]=min(abs(t2-output2));
[~,q3]=min(abs(t3-output2));
% t1=output2:output2+3;
nn1=find(t1==t1(q1,1));%n1是P波到时的后一个点
nn2=find(t2==t2(q2,1));%n1是P波到时的后一个点
nn3=find(t3==t3(q3,1));%n1是P波到时的后一个点



%%%峰值加速度的波形%%%%%
k=t4/dt;
Xg(1)=0; X(1)=0; 
for i=2:1:n1-1
    Xg(i)=0.999*Xg(i-1)+(inputA1(i-1)+inputA1(i+1)).^2;
    X(i)=0.999*X(i-1)+(2*inputA1(i)).^2;
    r=sqrt(Xg(i)/X(i));
    w(i)=real((1/(2*pi*dt))*acos(r));
end
%%%峰值加速度的波形%%%%%
omg_ud=w(nn1+1:nn1+k);

end

