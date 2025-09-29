clear;clc;
cd_UD = '\Data\';%�����ڶ����ļ�����
cd_EW = '\Data\EW\';%������һ���ļ�����
cd_NS = '\Data\NS\';%������һ���ļ�����

cd_move_UD = '\Data\pick\UD\';
cd_move_EW = '\Data\pick\EW\';
cd_move_NS = '\Data\pick\NS\';
mkdir(cd_move_UD);
mkdir(cd_move_EW);
mkdir(cd_move_NS);


g=dir(cd_UD);
w=g(3:end);                   %��ȡ�ļ���
[n1,~]=size(w);              %�ļ��и���

 %xiang guan can shu
 dat_info{1,1}='��ֵλ��Pd';dat_info{1,2}='��ֵ�ٶ�Pv';dat_info{1,3}='��ֵ���ٶ�Pa';
 dat_info{1,4}='�ٶ�ƽ������IV2';dat_info{1,5}='cav������ud���ۻ�';dat_info{1,6}='caa������ud���ۻ�';
 dat_info{1,7}='cad������ud���ۻ�';dat_info{1,8}='CAV�ۻ������ٶ�(������)';
 dat_info{1,9}='�ۻ������仯��DE';dat_info{1,10}='�������tc*Pd';dat_info{1,11}='��������tc';
 dat_info{1,12}='��ֵ��Tva'; dat_info{1,13}='��';dat_info{1,14}='���о�km';dat_info{1,15}='��Դ��km';
 dat_info{1,16}='��γ';dat_info{1,17}='��';dat_info{1,18}='��Դ���';dat_info{1,19}='�����';
 
 dat_info{1,20}='�����¼����';dat_info{1,21}='̨γ'; dat_info{1,22}='̨��';dat_info{1,23}='����Ƶ��';
 dat_info{1,24}='�����ٶ�';dat_info{1,25}='��ʱ';dat_info{1,26}='P����ʱ'; 
 

 for f1 = 1:1:n1               %��1���ļ�����ѭ�������¼�
    ReN=strcat(cd_UD,w(f1).name,'\','*.UD');    %Ѱ�Ҽ�¼������
    ReC=dir(ReN);           %ͳ��һ�����ٸ���¼
    [n2,~]=size(ReC);      %�Լ�¼��������
    
    k=2;
    for f2 = 1:1:n2           %��2���ļ�����ѭ�������¼
        Rec_UD=strcat(cd_UD,w(f1).name,'\',ReC(f2).name);
        str_UD = Rec_UD;
        %----%%%%read data %%%%--
        [DataUD,nt,t,ela,elo,sla,slo,M,sf,depth,Dir,Sta_H,Re,Maxacc,HD] = read_acc_data(str_UD);
        %Re���о�
        dt=1/nt; %ʱ������ntΪ����Ƶ��
        
        P_Acc = max(abs(DataUD));                         %��ֵ���ٶ�
        t_all=max(t);                                                %��¼ʱ��
        
        Rec_EW=strcat(cd_EW,ReC(f2).name(1:end-3),".EW");
        str_EW=Rec_EW;
        [DataEW,nt,t] = read_acc_data_EW_NS(str_EW);
        
        Rec_NS=strcat(cd_NS,ReC(f2).name(1:end-3),".NS");
        str_NS=Rec_NS;
        [DataNS,nt,t] = read_acc_data_EW_NS(str_NS);
        
        %�˹��޸�P����ʱ
        [a1,a2] = textread('\P-wave arrival.txt','%s%f','headerlines',0);%a1�ļ�����a2 P����ʱ
        p_num=size(a1,1);
        for p_n_1=1:p_num
            if  strcmp(strcat(a1{p_n_1,1},'.UD'),ReC(f2).name)==1 %,'.UD'
                P_a_t_r=round(a2(p_n_1)*nt)/nt;%p����ʱ 
                [SNR] = X_Z_B(DataUD,P_a_t_r,dt);
                [tc,Pd,TP,cad]= xzb_d_cs(DataUD,P_a_t_r,dt);
           
                [Pv,Pa,IV2,Tva,cav,caa,DE]= xzb_x_cs(DataUD,P_a_t_r,dt);
                [a_ud,t2]= bate_lvbo(DataUD,P_a_t_r,dt);%�˲���ʱ�䴰�ڵļ��ٶ�,t2��ʱ�䴰,UDfangxiang
                [a_ew,t2]= bate_lvbo(DataEW,P_a_t_r,dt);%�˲���ʱ�䴰�ڵļ��ٶ�,EWfangxiang
                [a_ns,t2]= bate_lvbo(DataNS,P_a_t_r,dt);%�˲���ʱ�䴰�ڵļ��ٶ�,NSfangxiang
                CAA= sqrt(((a_ud).^2)+((a_ew).^2)+((a_ns).^2));%ʱ�䴰��������ϳɼ��ٶ�
                CAV=max(cumtrapz(t2,abs(CAA)));%�ۻ������ٶȣ��������
                
            end
        end
        

        movefile(Rec_UD,cd_move_UD);
        movefile(Rec_EW,cd_move_EW);
        movefile(Rec_NS,cd_move_NS);
        
        
        dat_info{k,1}=Pd;dat_info{k,2}=Pv;dat_info{k,3}=Pa;
        dat_info{k,4}=IV2;dat_info{k,5}=cav;dat_info{k,6}=caa;
        dat_info{k,7}=cad;dat_info{k,8}=CAV;
        dat_info{k,9}=DE;dat_info{k,10}=TP;dat_info{k,11}=tc;
        dat_info{k,12}=Tva;dat_info{k,13}=M;dat_info{k,14}=Re;dat_info{k,15}=HD; 
        dat_info{k,16}=ela;dat_info{k,17}=elo;dat_info{k,18}=depth;dat_info{k,19}=SNR;
        
        dat_info{k,20}=ReC(f2).name;dat_info{k,21}=sla;dat_info{k,22}=slo;dat_info{k,23}=nt;
        dat_info{k,24}=P_Acc;dat_info{k,25}=t_all;dat_info{k,26}=P_a_t_r;
               
        
        k=k+1;
        %pause
    end
    
end
xlswrite('\SVM-Mag features.xlsx',dat_info);
