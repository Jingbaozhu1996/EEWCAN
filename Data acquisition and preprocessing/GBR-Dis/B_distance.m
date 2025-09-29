%�����о����B,y=Bt*exp(-At);
%������˹�˲�
function [B]= B_distance(inputA,output2,dt,T)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%inputA ������ٶ�
%output2 P����ʱ
%dt ʱ����
[m n]=size(inputA);%������ٶȵľ���ά��
for j=1:n
    t(j,1)=(j-1)*dt;
end
t4=T;%��ʾʱ�䴰����s
v=cumtrapz(t,inputA);%�ٶȣ����ǶԼ��ٶ������
d=cumtrapz(t,v);%λ�ƣ����Ƕ��ٶ������
%----------------4��0.075HZ-3HZ��ͨ�˲�-----------------------------------------%
                                                      %�˲��׶�
fmax=45;      fmin=40;
cutup=fmax/(1/(2*dt));
cutdown=fmin/(1/(2*dt));
wn=[cutdown, cutup]; 
[b,a]=butter(2,wn);                                               %4���˲�
%----------------4��0.075HZ-3HZ��ͨ�˲�-----------------------------------------%


%----------------4��0.075HZ��ͨ�˲�-----------------------------------------%
% Fs = 1/dt;
% fl= 0.075; % low cutoff frequency for highpass filter, depending on you
% ftype = 'high';
% f_order = 4;% order
% [b,a] = butter(f_order,fl/(Fs/2),ftype);% f_order order butterworth highpass filter     
%----------------4��0.075HZ��ͨ�˲�-----------------------------------------%
inputA=filter(b,a,inputA);%�˲���ļ��ٶ�
% d=filter(b,a,d);%�˲����λ��
% v=filter(b,a,v);%�˲�����ٶ�

y=abs(inputA)+0.001; %%y=|ACC|+e

[~,q]=min(abs(t-output2));
% t1=output2:output2+3;
n1=find(t==t(q,1));%n1��P����ʱ�ĵ�
% v1=v(n1:);
% d1=d(n1:);
% a1=inputA1(n1:n1+t4/0.01);
% t2=0:dt:t4;
% v1=cumtrapz(t2,a1);
% d1=cumtrapz(t2,v1);
%n1=n1-2;
% v1=v(n1:n1+t4/dt);%ʱ�䴰�ڵ��ٶ�
% d1=d(n1:n1+t4/dt);%ʱ�䴰�ڵ�λ��
% a1=y(n1:n1+t4/dt);%ʱ�䴰�ڵļ��ٶ�
MM=[];
YY=[];
tt=0.01;
ttt=[];
a=0;

for i=1:t4/dt
    MM(i)=max(y(n1:n1+i-1));  %%����������ϵĵ�
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

%%%���ư���ͼ
% AA=b(2);A=-AA;
% 
% subplot(2,1,1);
% [p_t,p_d] = getp_t_p_d(output2,1/dt,max(t),inputA,t,1,4);
% semilogy(p_t ,abs(p_d));
% axis([min(p_t),max(p_t),-inf,1000]);
% hold on
% P_Aplot=[output2,max(abs(p_d)); output2,0.000001];   %��P����ʱ��
% semilogy(P_Aplot(:,1),P_Aplot(:,2),'color','r','linewidth',0.1); %
% 
% hold on
% semilogy(ttt+output2-0.01,MM);
% 
% hold on
% yb=ttt.*exp(b(1)-A.*ttt);
% semilogy(ttt+output2-0.01,yb);
% title(name_data(1:end-3));  ylabel('���ٶ�(CM/sec/sec)');
% hold off
% 
% subplot(2,1,2);
% [p_t,p_d] = getp_t_p_d(output2,1/dt,max(t),inputA,t,1,4);
% plot(p_t ,abs(p_d));
% axis([min(p_t),max(p_t),0,max(abs(p_d))]);
% hold on
% P_Aplot=[output2,max(abs(p_d)); output2,0.000001];   %��P����ʱ��
% plot(P_Aplot(:,1),P_Aplot(:,2),'color','r','linewidth',0.1); %
% hold on
% plot(ttt+output2-0.01,MM);
% hold on
% yb=ttt.*exp(b(1)-A.*ttt);
% plot(ttt+output2-0.01,yb);
% hold off
end

