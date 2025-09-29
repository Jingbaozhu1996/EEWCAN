%����׿Խ���ڣ���ֵλ�ƣ���ֵ�ٶȣ���ֵ���ٶ�
%������˹�˲�
function [IV2,Tva,cav,caa,DE]= IV2_Tva_cav_caa_DE_lvbo_wave(inputA,inputB,inputC,output2,dt)
%UNTITLED Summary of this function goes here  
%   Detailed explanation goes here
%inputA ������ٶ�
%output2 P����ʱ
%dt ʱ����
[m n1]=size(inputA);%������ٶȵľ���ά��
for j=1:n1
    t1(j,1)=(j-1)*dt;
end

[m n2]=size(inputB);%������ٶȵľ���ά��
for j=1:n2
    t2(j,1)=(j-1)*dt;
end

[m n3]=size(inputC);%������ٶȵľ���ά��
for j=1:n3
    t3(j,1)=(j-1)*dt;
end


t4=3;%��ʾʱ�䴰����s
vA=cumtrapz(t1,inputA);%�ٶȣ����ǶԼ��ٶ������
dA=cumtrapz(t1,vA);%λ�ƣ����Ƕ��ٶ������

vB=cumtrapz(t2,inputB);%�ٶȣ����ǶԼ��ٶ������
dB=cumtrapz(t2,vB);%λ�ƣ����Ƕ��ٶ������

vC=cumtrapz(t3,inputC);%�ٶȣ����ǶԼ��ٶ������
dC=cumtrapz(t3,vC);%λ�ƣ����Ƕ��ٶ������
%----------------4��0.075HZ-3HZ��ͨ�˲�-----------------------------------------%
                                                      %�˲��׶�
fmax=3;      fmin=0.075;
cutup=fmax/(1/(2*dt));
cutdown=fmin/(1/(2*dt));
wn=[cutdown, cutup]; 
[b,a]=butter(4,wn);                                               %4���˲�
%----------------4��0.075HZ-3HZ��ͨ�˲�-----------------------------------------%


% %----------------4��0.075HZ��ͨ�˲�-----------------------------------------%
% Fs = 1/dt;
% fl= 0.075; % low cutoff frequency for highpass filter, depending on you
% ftype = 'high';
% f_order = 4;% order
% [b,a] = butter(f_order,fl/(Fs/2),ftype);% f_order order butterworth highpass filter     
% %----------------4��0.075HZ��ͨ�˲�-----------------------------------------%

inputA1=filter(b,a,inputA);%�˲���ļ��ٶ�
inputB1=filter(b,a,inputB);%�˲���ļ��ٶ�
inputC1=filter(b,a,inputC);%�˲���ļ��ٶ�

vA=filter(b,a,vA);%�˲�����ٶ�
vB=filter(b,a,vB);%�˲�����ٶ�
vC=filter(b,a,vC);%�˲�����ٶ�

dA=filter(b,a,dA);%�˲����λ��
dB=filter(b,a,dB);%�˲����λ��
dC=filter(b,a,dC);%�˲����λ��


[~,q1]=min(abs(t1-output2));
[~,q2]=min(abs(t2-output2));
[~,q3]=min(abs(t3-output2));
% t1=output2:output2+3;
nn1=find(t1==t1(q1,1));%n1��P����ʱ�ĺ�һ����
nn2=find(t2==t2(q2,1));%n1��P����ʱ�ĺ�һ����
nn3=find(t3==t3(q3,1));%n1��P����ʱ�ĺ�һ����

t2=0:dt:t4;

% v1=v(n1:n1+t4/dt);%ʱ�䴰�ڵ��ٶ�
% d1=d(n1:n1+t4/dt);%ʱ�䴰�ڵ�λ��
% a1=inputA1(n1:n1+t4/dt);%ʱ�䴰�ڵļ��ٶ�
% Pv=max(abs(v1));%P����ʱ��ʱ�䴰�ڵķ�ֵ�ٶ�




%%%%%%%%%%%%%%%%%%%%%%%
% t3=output2:dt:(output2+t4);
% IV2=trapz(t3,v1.^2);%�ٶ�ƽ������
%%%%%%%%%%%%%%%%%%%%%%%

%%%�ٶ�ƽ������IV2�Ĳ���%%%%%
k=t4/dt;
for i=1:1:k
    t3=output2:dt:(output2+dt*i);
    vud=vA(nn1:nn1+i);%ʱ�䴰�ڵ��ٶ�
    vew=vB(nn2:nn2+i);%ʱ�䴰�ڵ��ٶ�
    vns=vC(nn3:nn3+i);%ʱ�䴰�ڵ��ٶ�
    V_hecheng= sqrt(((vud).^2)+((vew).^2)+((vns).^2));%ʱ�䴰��������ϳɼ��ٶ�
    IV2(i)=trapz(t3,V_hecheng.^2);%�ٶ�ƽ������
end
%%%�ٶ�ƽ�����ֵĲ���%%%%%   

%%%��ֵ��Tva%%%%%%%%%%%%%%%%%%%%% 
k=t4/dt;
for i=1:1:k
    vud=abs(vA(nn1:nn1+i));
    vew=abs(vB(nn2:nn2+i));
    vns=abs(vC(nn3:nn3+i));
    Pv_hecheng= max(sqrt(((vud).^2)+((vew).^2)+((vns).^2)));%ʱ�䴰��������ϳɼ��ٶ�
    
    aud=abs(inputA1(nn1:nn1+i));
    aew=abs(inputB1(nn2:nn2+i));
    ans=abs(inputC1(nn3:nn3+i));
    Pa_hecheng= max(sqrt(((aud).^2)+((aew).^2)+((ans).^2)));%ʱ�䴰��������ϳɼ��ٶ�
    Tva(i)=2*pi*(Pv_hecheng/Pa_hecheng);%P����ʱ��ʱ�䴰�ڵķ�ֵ��
end
%%%��ֵ��Tva%%%%%%%%%%%%%%%%%%%%%  

%%%cav,caa%%%%%%%%%%%%%%%%%%%%%
k=t4/dt;
for i=1:1:k
    cavA=abs(vA(nn1:nn1+i));%cav�ۻ������ٶ�UD
    cavB=abs(vB(nn2:nn2+i));%cav�ۻ������ٶ�UD
    cavC=abs(vC(nn3:nn3+i));%cav�ۻ������ٶ�UD
    cav(i)= sum(sqrt(((cavA).^2)+((cavB).^2)+((cavC).^2)));%ʱ�䴰��������ϳɼ��ٶ�
    
    caaA=abs(inputA1(nn1:nn1+i));%caa�ۻ����Լ��ٶ�UD
    caaB=abs(inputB1(nn2:nn2+i));%caa�ۻ����Լ��ٶ�UD
    caaC=abs(inputC1(nn3:nn3+i));%caa�ۻ����Լ��ٶ�UD
    caa(i)= sum(sqrt(((caaA).^2)+((caaB).^2)+((caaC).^2)));%ʱ�䴰��������ϳɼ��ٶ�
end
%%%cav,caa%%%%%%%%%%%%%%%%%%%%%

%%%DE%%%%%%%%%%%%%%%%%%%%%
k=t4/dt;
for i=1:1:k
    vud=vA(nn1:nn1+i);
    vew=vB(nn2:nn2+i);
    vns=vC(nn3:nn3+i);
    v_hecheng= sqrt(((vud).^2)+((vew).^2)+((vns).^2));%ʱ�䴰��������ϳ��ٶ�
    
    aud=inputA1(nn1:nn1+i); 
    aew=inputB1(nn2:nn2+i); 
    ans=inputC1(nn3:nn3+i);
    a_heceng= sqrt(((aud).^2)+((aew).^2)+((ans).^2));%ʱ�䴰��������ϳɼ��ٶ�

    DI=log(abs(a_heceng.*v_hecheng));
    DE(i)=max(DI);%�ۻ������仯��
end
%%%DE%%%%%%%%%%%%%%%%%%%%%



end

