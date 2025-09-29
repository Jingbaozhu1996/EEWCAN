%计算卓越周期，峰值位移，峰值速度，峰值加速度
%巴特沃斯滤波
function [Tvaud]= Tvaud_lvbo_wave(inputA,output2,dt)
%UNTITLED Summary of this function goes here  
%   Detailed explanation goes here
%inputA 输入加速度
%output2 P波到时
%dt 时间间隔
[m n1]=size(inputA);%输入加速度的矩阵维度
for j=1:n1
    t1(j,1)=(j-1)*dt;
end



t4=3;%表示时间窗长度s
vA=cumtrapz(t1,inputA);%速度，这是对加速度求积分
dA=cumtrapz(t1,vA);%位移，这是对速度求积分


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

vA=filter(b,a,vA);%滤波后的速度

dA=filter(b,a,dA);%滤波后的位移


[~,q1]=min(abs(t1-output2));

% t1=output2:output2+3;
nn1=find(t1==t1(q1,1));%n1是P波到时的后一个点


t2=0:dt:t4;



%%%峰值比Tva%%%%%%%%%%%%%%%%%%%%% 
k=t4/dt;
for i=1:1:k
    vud=abs(vA(nn1:nn1+i));
    
    aud=abs(inputA1(nn1:nn1+i));

    Tvaud(i)=2*pi*(vud/aud);%P波到时后，时间窗内的峰值比
end
%%%峰值比Tva%%%%%%%%%%%%%%%%%%%%%  


end

