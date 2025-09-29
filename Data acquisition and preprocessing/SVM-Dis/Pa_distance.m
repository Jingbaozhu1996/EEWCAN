%计算卓越周期，峰值位移，峰值速度，峰值加速度
%巴特沃斯滤波
function [Pa]= Pa_distance(inputA,output2,dt,T,J)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%inputA 输入加速度
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
inputA1=filter(b,a,inputA);%滤波后的加速度
d=filter(b,a,d);%滤波后的位移
v=filter(b,a,v);%滤波后的速度
[~,q]=min(abs(t-output2));
% t1=output2:output2+3;
n1=find(t==t(q,1));%n1是P波到时的后一个点
% v1=v(n1:);
% d1=d(n1:);
% a1=inputA1(n1:n1+t4/0.01);
t2=0:dt:t4;
% v1=cumtrapz(t2,a1);
% d1=cumtrapz(t2,v1);
v1=v(n1:n1+t4/dt);%时间窗内的速度
d1=d(n1:n1+t4/dt);%时间窗内的位移
a1=inputA1(n1:n1+t4/dt);%时间窗内的加速度
Pa=max(abs(a1));%P波到时后，时间窗内的峰值位移


end

