%计算卓越周期，峰值位移，峰值速度，峰值加速度
%巴特沃斯滤波
function [tcud,TPud]= tcud_TPud_lvbo_wave(inputA,output2,dt)
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
% fmax=3;      fmin=0.075;
% cutup=fmax/(1/(2*dt));
% cutdown=fmin/(1/(2*dt));
% wn=[cutdown, cutup]; 
% [b,a]=butter(4,wn);                                               %4阶滤波
%----------------4阶0.075HZ-3HZ带通滤波-----------------------------------------%


%----------------4阶0.075HZ高通滤波-----------------------------------------%
Fs = 1/dt;
fl= 0.075; % low cutoff frequency for highpass filter, depending on you
ftype = 'high';  %low
f_order = 4;% order
[b,a] = butter(f_order,fl/(Fs/2),ftype);% f_order order butterworth highpass filter     
%----------------4阶0.075HZ高通滤波-----------------------------------------%
inputA1=filter(b,a,inputA);%滤波后的加速度


vA=filter(b,a,vA);%滤波后的速度


dA=filter(b,a,dA);%滤波后的位移


[~,q1]=min(abs(t1-output2));

% t1=output2:output2+3;
nn1=find(t1==t1(q1,1));%n1是P波到时的后一个点

% v1=v(n1:);
% d1=d(n1:);
% a1=inputA1(n1:n1+t4/0.01);
t2=0:dt:t4;
% v1=cumtrapz(t2,a1);
% d1=cumtrapz(t2,v1);
%    v1=v(n1:n1+t4/dt);%时间窗内的速度
%    d1=d(n1:n1+t4/dt);%时间窗内的位移
%    a1=inputA1(n1:n1+t4/dt);%时间窗内的加速度
%    Pd=max(abs(d1));%P波到时后，时间窗内的峰值位移
% % r1=quad(inline('(v.^2)'),0,4);
% % r2=quad(inline('(d.^2)'),0,4);
%    r1=trapz(t2,v1.^2);
%    r2=trapz(t2,d1.^2);
%    r=r1/r2;
%tc=2*pi/sqrt(r);%特征周期
%TP=tc*Pd;%构造参数


%------------计算tc三分向参数波形---------%
k=t4/dt;
for i=1:1:k   
    t2=output2:dt:(output2+dt*i);
    vud=vA(nn1:nn1+i);%时间窗内的速度

    
%    v1=v(n1:n1+i);%时间窗内的速度
    dud=dA(nn1:nn1+i);%时间窗内的位移

%    d1=d(n1:n1+i);%时间窗内的位移
    r1=trapz(t2,vud.^2);
    r2=trapz(t2,dud.^2);
    r=r1/r2;
    tcud(i)=2*pi/sqrt(r);%特征周期
end
%------------计算tc参数波形---------%




%------------计算TP参数波形---------%
k=t4/dt;
for i=1:1:k   
    t2=output2:dt:(output2+dt*i);
    vud=vA(nn1:nn1+i);%时间窗内的速度

    
%    v1=v(n1:n1+i);%时间窗内的速度
    dud=dA(nn1:nn1+i);%时间窗内的位移

%    d1=d(n1:n1+i);%时间窗内的位移
    r1=trapz(t2,vud.^2);
    r2=trapz(t2,dud.^2);
    r=r1/r2;
    tcc=2*pi/sqrt(r);%特征周期
    

    Pd= max(abs(dud));%时间窗内三分向合成位移   
    
    TPud(i)=tcc*Pd %构造参数
end
%------------计算TP参数波形---------%



end

