%求震中距参数B,y=Bt*exp(-At);
%巴特沃斯滤波
function [B]= B_distance(inputA,output2,dt,T)
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
[b,a]=butter(2,wn);                                               %4阶滤波
%----------------4阶0.075HZ-3HZ带通滤波-----------------------------------------%


%----------------4阶0.075HZ高通滤波-----------------------------------------%
% Fs = 1/dt;
% fl= 0.075; % low cutoff frequency for highpass filter, depending on you
% ftype = 'high';
% f_order = 4;% order
% [b,a] = butter(f_order,fl/(Fs/2),ftype);% f_order order butterworth highpass filter     
%----------------4阶0.075HZ高通滤波-----------------------------------------%
inputA=filter(b,a,inputA);%滤波后的加速度
% d=filter(b,a,d);%滤波后的位移
% v=filter(b,a,v);%滤波后的速度

y=abs(inputA)+0.001; %%y=|ACC|+e

[~,q]=min(abs(t-output2));
% t1=output2:output2+3;
n1=find(t==t(q,1));%n1是P波到时的点
% v1=v(n1:);
% d1=d(n1:);
% a1=inputA1(n1:n1+t4/0.01);
% t2=0:dt:t4;
% v1=cumtrapz(t2,a1);
% d1=cumtrapz(t2,v1);
%n1=n1-2;
% v1=v(n1:n1+t4/dt);%时间窗内的速度
% d1=d(n1:n1+t4/dt);%时间窗内的位移
% a1=y(n1:n1+t4/dt);%时间窗内的加速度
MM=[];
YY=[];
tt=0.01;
ttt=[];
a=0;

for i=1:t4/dt
    MM(i)=max(y(n1:n1+i-1));  %%计算包络线上的点
    YY(i)=log(MM(i)/tt);    
    tt=tt+0.01;
    ttt(i)=a+0.01;
    a=a+0.01;
end

ttt=ttt';
YY=YY';
X=[ones(length(ttt),1),ttt];Y=YY; %y=lnB+At
[b,bint,r,rint,stats]=regress(Y,X); 
B=exp(b(1));

%%%绘制包络图
% AA=b(2);A=-AA;
% 
% subplot(2,1,1);
% [p_t,p_d] = getp_t_p_d(output2,1/dt,max(t),inputA,t,1,4);
% semilogy(p_t ,abs(p_d));
% axis([min(p_t),max(p_t),-inf,1000]);
% hold on
% P_Aplot=[output2,max(abs(p_d)); output2,0.000001];   %画P波到时线
% semilogy(P_Aplot(:,1),P_Aplot(:,2),'color','r','linewidth',0.1); %
% 
% hold on
% semilogy(ttt+output2-0.01,MM);
% 
% hold on
% yb=ttt.*exp(b(1)-A.*ttt);
% semilogy(ttt+output2-0.01,yb);
% title(name_data(1:end-3));  ylabel('加速度(CM/sec/sec)');
% hold off
% 
% subplot(2,1,2);
% [p_t,p_d] = getp_t_p_d(output2,1/dt,max(t),inputA,t,1,4);
% plot(p_t ,abs(p_d));
% axis([min(p_t),max(p_t),0,max(abs(p_d))]);
% hold on
% P_Aplot=[output2,max(abs(p_d)); output2,0.000001];   %画P波到时线
% plot(P_Aplot(:,1),P_Aplot(:,2),'color','r','linewidth',0.1); %
% hold on
% plot(ttt+output2-0.01,MM);
% hold on
% yb=ttt.*exp(b(1)-A.*ttt);
% plot(ttt+output2-0.01,yb);
% hold off
end

