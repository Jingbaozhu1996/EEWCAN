%计算卓越周期，峰值位移，峰值速度，峰值加速度
%巴特沃斯滤波
function [IV2,Tva,cav,caa,DE]= IV2_Tva_cav_caa_DE_lvbo_wave(inputA,inputB,inputC,output2,dt)
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


% %----------------4阶0.075HZ高通滤波-----------------------------------------%
% Fs = 1/dt;
% fl= 0.075; % low cutoff frequency for highpass filter, depending on you
% ftype = 'high';
% f_order = 4;% order
% [b,a] = butter(f_order,fl/(Fs/2),ftype);% f_order order butterworth highpass filter     
% %----------------4阶0.075HZ高通滤波-----------------------------------------%

inputA1=filter(b,a,inputA);%滤波后的加速度
inputB1=filter(b,a,inputB);%滤波后的加速度
inputC1=filter(b,a,inputC);%滤波后的加速度

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

% v1=v(n1:n1+t4/dt);%时间窗内的速度
% d1=d(n1:n1+t4/dt);%时间窗内的位移
% a1=inputA1(n1:n1+t4/dt);%时间窗内的加速度
% Pv=max(abs(v1));%P波到时后，时间窗内的峰值速度




%%%%%%%%%%%%%%%%%%%%%%%
% t3=output2:dt:(output2+t4);
% IV2=trapz(t3,v1.^2);%速度平方积分
%%%%%%%%%%%%%%%%%%%%%%%

%%%速度平方积分IV2的波形%%%%%
k=t4/dt;
for i=1:1:k
    t3=output2:dt:(output2+dt*i);
    vud=vA(nn1:nn1+i);%时间窗内的速度
    vew=vB(nn2:nn2+i);%时间窗内的速度
    vns=vC(nn3:nn3+i);%时间窗内的速度
    V_hecheng= sqrt(((vud).^2)+((vew).^2)+((vns).^2));%时间窗内三分向合成加速度
    IV2(i)=trapz(t3,V_hecheng.^2);%速度平方积分
end
%%%速度平方积分的波形%%%%%   

%%%峰值比Tva%%%%%%%%%%%%%%%%%%%%% 
k=t4/dt;
for i=1:1:k
    vud=abs(vA(nn1:nn1+i));
    vew=abs(vB(nn2:nn2+i));
    vns=abs(vC(nn3:nn3+i));
    Pv_hecheng= max(sqrt(((vud).^2)+((vew).^2)+((vns).^2)));%时间窗内三分向合成加速度
    
    aud=abs(inputA1(nn1:nn1+i));
    aew=abs(inputB1(nn2:nn2+i));
    ans=abs(inputC1(nn3:nn3+i));
    Pa_hecheng= max(sqrt(((aud).^2)+((aew).^2)+((ans).^2)));%时间窗内三分向合成加速度
    Tva(i)=2*pi*(Pv_hecheng/Pa_hecheng);%P波到时后，时间窗内的峰值比
end
%%%峰值比Tva%%%%%%%%%%%%%%%%%%%%%  

%%%cav,caa%%%%%%%%%%%%%%%%%%%%%
k=t4/dt;
for i=1:1:k
    cavA=abs(vA(nn1:nn1+i));%cav累积绝对速度UD
    cavB=abs(vB(nn2:nn2+i));%cav累积绝对速度UD
    cavC=abs(vC(nn3:nn3+i));%cav累积绝对速度UD
    cav(i)= sum(sqrt(((cavA).^2)+((cavB).^2)+((cavC).^2)));%时间窗内三分向合成加速度
    
    caaA=abs(inputA1(nn1:nn1+i));%caa累积绝对加速度UD
    caaB=abs(inputB1(nn2:nn2+i));%caa累积绝对加速度UD
    caaC=abs(inputC1(nn3:nn3+i));%caa累积绝对加速度UD
    caa(i)= sum(sqrt(((caaA).^2)+((caaB).^2)+((caaC).^2)));%时间窗内三分向合成加速度
end
%%%cav,caa%%%%%%%%%%%%%%%%%%%%%

%%%DE%%%%%%%%%%%%%%%%%%%%%
k=t4/dt;
for i=1:1:k
    vud=vA(nn1:nn1+i);
    vew=vB(nn2:nn2+i);
    vns=vC(nn3:nn3+i);
    v_hecheng= sqrt(((vud).^2)+((vew).^2)+((vns).^2));%时间窗内三分向合成速度
    
    aud=inputA1(nn1:nn1+i); 
    aew=inputB1(nn2:nn2+i); 
    ans=inputC1(nn3:nn3+i);
    a_heceng= sqrt(((aud).^2)+((aew).^2)+((ans).^2));%时间窗内三分向合成加速度

    DI=log(abs(a_heceng.*v_hecheng));
    DE(i)=max(DI);%累积能量变化率
end
%%%DE%%%%%%%%%%%%%%%%%%%%%



end

