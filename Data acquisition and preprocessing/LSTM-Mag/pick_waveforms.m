clear;clc;
cd_UD ='\Data\';%�����ڶ����ļ�����
cd_EW ='\Data\EW\';%������һ���ļ�����
cd_NS ='\Data\NS\';%������һ���ļ�����

cd_move_UD = '\Data\pick\UD\';
cd_move_EW = '\Data\pick\EW\';
cd_move_NS = '\Data\pick\NS\';
mkdir(cd_move_UD);
mkdir(cd_move_EW);
mkdir(cd_move_NS);


g=dir(cd_UD);
w=g(3:end);                   %��ȡ�ļ���
[n1,~]=size(w);              %�ļ��и���
 dat_info{1,1}='�����¼����';dat_info{1,2}='��γ';dat_info{1,3}='��';dat_info{1,4}='��Դ���'; 
 dat_info{1,5}='̨γ'; dat_info{1,6}='̨��';dat_info{1,7}='����Ƶ��';
 dat_info{1,8}='�����ٶ�';dat_info{1,9}='��ʱ';dat_info{1,10}='��Դ��km';
 %xiang guan can shu
 dat_info{1,11}='��';dat_info{1,12}='���о�km';
 dat_info{1,13}='PGA_MMI';dat_info{1,14}='PGV_MMI';
 dat_info{1,15}='PGA_TIME';dat_info{1,16}='PGV_TIME';
 dat_info{1,17}='�����';dat_info{1,18}='P����ʱ';
 
 
 ACC_UD=zeros(10856,300);ACC_EW=zeros(10856,300);ACC_NS=zeros(10856,300);
 VEL_UD=zeros(10856,300);VEL_EW=zeros(10856,300);VEL_NS=zeros(10856,300);
 DIS_UD=zeros(10856,300);DIS_EW=zeros(10856,300);DIS_NS=zeros(10856,300);
 
  
 z=1;
 k=2;
 for f1 = 1:1:n1               %��1���ļ�����ѭ�������¼�
    ReN=strcat(cd_UD,w(f1).name,'\','*.UD');    %Ѱ�Ҽ�¼������
    ReC=dir(ReN);           %ͳ��һ�����ٸ���¼
    [n2,~]=size(ReC);      %�Լ�¼��������
    
    for f2 = 1:1:n2           %��2���ļ�����ѭ�������¼
        Rec_UD=strcat(cd_UD,w(f1).name,'\',ReC(f2).name);
        str_UD = Rec_UD;
        %----%%%%read data%%%%--
        [DataUD,nt,t,ela,elo,sla,slo,M,sf,depth,Dir,Sta_H,Re,Maxacc,HD] = read_acc_data(str_UD);
        %Re���о�
        dt=1/nt;%ʱ������ntΪ����Ƶ��
        
        
        t_all=max(t);                             %��¼ʱ��
        P_Acc = max(abs(DataUD)); 
        
        Rec_EW=strcat(cd_EW,ReC(f2).name(1:end-3),".EW");
        str_EW=Rec_EW;
        [DataEW,~,~] = read_acc_data_EW_NS(str_EW);
        %[Gaotie_Acc_EW,Gaotie_Vel_EW]= Gaotie_Jisuan_PGA_PGV(DataEW,dt); %���� EW����PGA��PGV
        [MMI_ACC_EW,MMI_VEL_EW] = MMI_Jisuan_ACC_VEL(DataEW,dt);  %MMI EW����PGA��PGV
        
        
        Rec_NS=strcat(cd_NS,ReC(f2).name(1:end-3),".NS");
        str_NS=Rec_NS;
        [DataNS,nt,t] = read_acc_data_EW_NS(str_NS);
        %[Gaotie_Acc_NS,Gaotie_Vel_NS] = Gaotie_Jisuan_PGA_PGV(DataNS,dt); %���� NS����PGA��PGV
        [MMI_ACC_NS,MMI_VEL_NS] = MMI_Jisuan_ACC_VEL(DataNS,dt);  %MMI NS����PGA��PGV
        
        %�ϳ�PGA
        PGA2_MMI = max(sqrt(((MMI_ACC_EW).^2)+((MMI_ACC_NS).^2)));%MMIʱ�䴰��ˮƽ����ϳɼ��ٶ�
        %�ϳ�PGV
        PGV2_MMI= max(sqrt(((MMI_VEL_EW).^2)+((MMI_VEL_NS).^2)));%MMIʱ�䴰��ˮƽ����ϳ��ٶ�
        
        %%����PGA��Ӧ��ʱ��      
        A_hecheng=sqrt(((MMI_ACC_EW).^2)+((MMI_ACC_NS).^2));
        pga_po=find(PGA2_MMI==A_hecheng);
        T_PGA=pga_po*dt;
        
        %%����PGV��Ӧ��ʱ��
        V_hecheng=sqrt(((MMI_VEL_EW).^2)+((MMI_VEL_NS).^2));
        pgv_po=find(PGV2_MMI==V_hecheng);
        T_PGV=pgv_po*dt;        
        

        %�˹��޸�P����ʱ
        [a1,a2] = textread('\P-wave arrival.txt','%s%f','headerlines',0);%a1�ļ�����a2 P����ʱ
        p_num=size(a1,1);
        for p_n_1=1:p_num
            if  strcmp(strcat(a1{p_n_1,1},'.UD'),ReC(f2).name)==1 %,'.UD'
                P_a_t_r=round(a2(p_n_1)*nt)/nt;%p����ʱ 
                [SNR]= X_Z_B(DataUD,P_a_t_r,dt);
                [d_ud,d_ew,d_ns]= raw_lvbo_dis(DataUD,DataEW,DataNS,P_a_t_r,dt,3);%0.075HZ��ͨ�˲�
                [v_ud,v_ew,v_ns]= raw_lvbo_vel(DataUD,DataEW,DataNS,P_a_t_r,dt,3);%0.075HZ��ͨ�˲�
                [a_ud,a_ew,a_ns]= raw_lvbo_acc(DataUD,DataEW,DataNS,P_a_t_r,dt,3);%0.075HZ��ͨ�˲�
               
                
            end
        end
        

        movefile(Rec_UD,cd_move_UD);
        movefile(Rec_EW,cd_move_EW);
        movefile(Rec_NS,cd_move_NS);
        
        
        
        ACC_UD(z,:)=a_ud;ACC_EW(z,:)=a_ew;ACC_NS(z,:)=a_ns;
        VEL_UD(z,:)=v_ud;VEL_EW(z,:)=v_ew;VEL_NS(z,:)=v_ns;
        DIS_UD(z,:)=d_ud;DIS_EW(z,:)=d_ew;DIS_NS(z,:)=d_ns;
        
        
        
        dat_info{k,1}=ReC(f2).name;dat_info{k,2}=ela;dat_info{k,3}=elo;dat_info{k,4}=depth; 
        dat_info{k,5}=sla;dat_info{k,6}=slo;dat_info{k,7}=nt;dat_info{k,8}=P_Acc;
        dat_info{k,9}=t_all;dat_info{k,10}=HD; 
        
        dat_info{k,11}=M;dat_info{k,12}=Re;
        dat_info{k,13}=PGA2_MMI;dat_info{k,14}=PGV2_MMI;
        dat_info{k,15}=T_PGA;dat_info{k,16}=T_PGV;
        dat_info{k,17}=SNR;dat_info{k,18}=P_a_t_r;
        
        k=k+1;
        z=z+1;
        %pause
    end
    
end
xlswrite('\LSTM-Mag\xinxi.xlsx',dat_info);

save('\LSTM-Mag\ACC_UD.mat','ACC_UD'); 
save('\LSTM-Mag\ACC_EW.mat','ACC_EW');
save('\LSTM-Mag\ACC_NS.mat','ACC_NS');

save('\LSTM-Mag\VEL_UD.mat','VEL_UD');
save('\LSTM-Mag\VEL_EW.mat','VEL_EW');
save('\LSTM-Mag\VEL_NS.mat','VEL_NS');

save('\LSTM-Mag\DIS_UD.mat','DIS_UD'); 
save('\LSTM-Mag\DIS_EW.mat','DIS_EW');
save('\LSTM-Mag\DIS_NS.mat','DIS_NS');

