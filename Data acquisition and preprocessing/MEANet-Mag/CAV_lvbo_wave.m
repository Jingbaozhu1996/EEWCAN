%巴特沃思滤波
function [CAV]= CAV_lvbo_wave(inputA,inputB,inputC,output2,dt)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%inputA 输入UD加速度
%inputB 输入EW加速度
%inputC 输入NS加速度
%output2 P波到时
%dt 时间间隔
%a1滤波后时间窗内的加速度
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


% %----------------4阶0.075HZ高通滤波-----------------------------------------%
% Fs = 1/dt;
% fl= 0.075; % low cutoff frequency for highpass filter, depending on you
% ftype = 'high';
% f_order = 4;% order
% [b,a] = butter(f_order,fl/(Fs/2),ftype);% f_order order butterworth highpass filter     
% %----------------4阶0.075HZ高通滤波-----------------------------------------%


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

% v1=v(n1:);
% d1=d(n1:);
% a1=inputA1(n1:n1+t4/0.01);
k=t4/dt;
for i=1:1:k
    t2=output2:dt:(output2+dt*i);
    aud=inputA1(nn1:nn1+i);%ud加速度
    aew=inputB1(nn2:nn2+i);%ew加速度
    ans=inputC1(nn3:nn3+i);%ew加速度
    CAA= sqrt(((aud).^2)+((aew).^2)+((ans).^2));%时间窗内三分向合成加速度
    CAV(i)=max(cumtrapz(t2,abs(CAA)));%累积绝对速度，三分向的
end


end