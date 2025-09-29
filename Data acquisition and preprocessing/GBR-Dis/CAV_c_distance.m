%求震中距参数D,y=Dt.根据CAV
%巴特沃斯滤波
function [CAV_c]= CAV_c_distance(inputA,inputB,inputC,output2,dt,T,J)
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
fmax=40;      fmin=35;
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

%CAA= sqrt(((a_ud).^2)+((a_ew).^2)+((a_ns).^2));%时间窗内三分向合成加速度
yy=cumtrapz(t,abs(syn_acc));%累积绝对速度，三分向的
%y=syn_acc+0.001;

% yy=[];
% yy(1)=y(1);
% for i=2:n
%     yy(i)=yy(i-1)+y(i);
% end

[~,q]=min(abs(t-output2));
% t1=output2:output2+3;
n1=find(t==t(q,1));%n1是P波到时的后一个点
% v1=v(n1:);
% d1=d(n1:);
% a1=inputA1(n1:n1+t4/0.01);
t2=0:dt:t4;
% v1=cumtrapz(t2,a1);
% d1=cumtrapz(t2,v1);
n1=n1-2;
%v1=v(n1:n1+t4/dt);%时间窗内的速度
%d1=d(n1:n1+t4/dt);%时间窗内的位移
%a1=y(n1:n1+t4/dt);%时间窗内的加速度
MM=yy(n1-1:n1-2+t4/dt);
% for i=1:t4/dt
%     MM(i)=yy(n1+1:n1+i);
% end
ttt=[];
a=0;
for k=1:t4/dt
    ttt(k)=a+0.01;
    a=a+0.01;   
end
ttt=ttt';
MM=MM';
X=[ones(length(ttt),1),ttt];Y=[MM];
[b,bint,r,rint,stats]=regress(Y,X); 
CAV_c=b(2);
end

