%计算MMI烈度PGA\PGV计算
%巴特沃斯滤波
function [MMI_ACC,MMI_VEL]= MMI_Jisuan_ACC_VEL(inputA,dt)
%UNTITLED Summary of this function goes here
%Detailed explanation goes here
%inputA 输入加速度
%output2 P波到时
%dt 时间间隔
[~, n]=size(inputA);%输入加速度的矩阵维度
for j=1:n
    t(j,1)=(j-1)*dt;
end

v=cumtrapz(t,inputA);%速度，这是对加速度求积分

%----------------4阶0.075HZ-3HZ带通滤波-----------------------------------------%
                                                      %滤波阶段
% fmax=0.1;      fmin=10;
% cutup=fmax/(1/(2*dt));
% cutdown=fmin/(1/(2*dt));
% wn=[cutdown, cutup]; 
% [b,a]=butter(4,wn);                                               %4阶滤波
%----------------4阶0.075HZ-3HZ带通滤波-----------------------------------------%


% %----------------4阶0.075HZ高通滤波-----------------------------------------%
Fs = 1/dt;
fl= 0.075; % low cutoff frequency for highpass filter, depending on you
ftype = 'high';
f_order = 4;% order
[b,a] = butter(f_order,fl/(Fs/2),ftype);% f_order order butterworth highpass filter     
% %----------------4阶0.075HZ高通滤波-----------------------------------------%

inputA1=filter(b,a,inputA);%滤波后的加速度
MMI_ACC=inputA1;%峰值加速度

v=filter(b,a,v);%滤波后的速度
MMI_VEL=v;%峰值速度



end

