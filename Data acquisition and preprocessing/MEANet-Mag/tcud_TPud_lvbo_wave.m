%����׿Խ���ڣ���ֵλ�ƣ���ֵ�ٶȣ���ֵ���ٶ�
%������˹�˲�
function [tcud,TPud]= tcud_TPud_lvbo_wave(inputA,output2,dt)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%inputA ������ٶ�
%output2 P����ʱ
%dt ʱ����
[m n1]=size(inputA);%������ٶȵľ���ά��
for j=1:n1
    t1(j,1)=(j-1)*dt;
end



t4=3;%��ʾʱ�䴰����s
vA=cumtrapz(t1,inputA);%�ٶȣ����ǶԼ��ٶ������
dA=cumtrapz(t1,vA);%λ�ƣ����Ƕ��ٶ������


%----------------4��0.075HZ-3HZ��ͨ�˲�-----------------------------------------%
                                                      %�˲��׶�
% fmax=3;      fmin=0.075;
% cutup=fmax/(1/(2*dt));
% cutdown=fmin/(1/(2*dt));
% wn=[cutdown, cutup]; 
% [b,a]=butter(4,wn);                                               %4���˲�
%----------------4��0.075HZ-3HZ��ͨ�˲�-----------------------------------------%


%----------------4��0.075HZ��ͨ�˲�-----------------------------------------%
Fs = 1/dt;
fl= 0.075; % low cutoff frequency for highpass filter, depending on you
ftype = 'high';  %low
f_order = 4;% order
[b,a] = butter(f_order,fl/(Fs/2),ftype);% f_order order butterworth highpass filter     
%----------------4��0.075HZ��ͨ�˲�-----------------------------------------%
inputA1=filter(b,a,inputA);%�˲���ļ��ٶ�


vA=filter(b,a,vA);%�˲�����ٶ�


dA=filter(b,a,dA);%�˲����λ��


[~,q1]=min(abs(t1-output2));

% t1=output2:output2+3;
nn1=find(t1==t1(q1,1));%n1��P����ʱ�ĺ�һ����

% v1=v(n1:);
% d1=d(n1:);
% a1=inputA1(n1:n1+t4/0.01);
t2=0:dt:t4;
% v1=cumtrapz(t2,a1);
% d1=cumtrapz(t2,v1);
%    v1=v(n1:n1+t4/dt);%ʱ�䴰�ڵ��ٶ�
%    d1=d(n1:n1+t4/dt);%ʱ�䴰�ڵ�λ��
%    a1=inputA1(n1:n1+t4/dt);%ʱ�䴰�ڵļ��ٶ�
%    Pd=max(abs(d1));%P����ʱ��ʱ�䴰�ڵķ�ֵλ��
% % r1=quad(inline('(v.^2)'),0,4);
% % r2=quad(inline('(d.^2)'),0,4);
%    r1=trapz(t2,v1.^2);
%    r2=trapz(t2,d1.^2);
%    r=r1/r2;
%tc=2*pi/sqrt(r);%��������
%TP=tc*Pd;%�������


%------------����tc�������������---------%
k=t4/dt;
for i=1:1:k   
    t2=output2:dt:(output2+dt*i);
    vud=vA(nn1:nn1+i);%ʱ�䴰�ڵ��ٶ�

    
%    v1=v(n1:n1+i);%ʱ�䴰�ڵ��ٶ�
    dud=dA(nn1:nn1+i);%ʱ�䴰�ڵ�λ��

%    d1=d(n1:n1+i);%ʱ�䴰�ڵ�λ��
    r1=trapz(t2,vud.^2);
    r2=trapz(t2,dud.^2);
    r=r1/r2;
    tcud(i)=2*pi/sqrt(r);%��������
end
%------------����tc��������---------%




%------------����TP��������---------%
k=t4/dt;
for i=1:1:k   
    t2=output2:dt:(output2+dt*i);
    vud=vA(nn1:nn1+i);%ʱ�䴰�ڵ��ٶ�

    
%    v1=v(n1:n1+i);%ʱ�䴰�ڵ��ٶ�
    dud=dA(nn1:nn1+i);%ʱ�䴰�ڵ�λ��

%    d1=d(n1:n1+i);%ʱ�䴰�ڵ�λ��
    r1=trapz(t2,vud.^2);
    r2=trapz(t2,dud.^2);
    r=r1/r2;
    tcc=2*pi/sqrt(r);%��������
    

    Pd= max(abs(dud));%ʱ�䴰��������ϳ�λ��   
    
    TPud(i)=tcc*Pd %�������
end
%------------����TP��������---------%



end

