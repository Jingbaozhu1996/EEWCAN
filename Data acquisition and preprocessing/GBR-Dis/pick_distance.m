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
w=g(3:end);                  %��ȡ�ļ���
[n1,~]=size(w);              %�ļ��и���

 dat_info{1,1}='�����';dat_info{1,2}='P����ʱ';dat_info{1,3}='M_type';
 dat_info{1,4}='�����¼����';dat_info{1,5}='��γ';dat_info{1,6}='��';dat_info{1,7}='��Դ���'; 
 dat_info{1,8}='��';dat_info{1,9}='̨γ'; dat_info{1,10}='̨��';dat_info{1,11}='����Ƶ��';
 dat_info{1,12}='���о�km';dat_info{1,13}='��Դ��km';
 
 %%%��ֵ�༰����������%%%%%
 dat_info{1,14}='Pdֵ';dat_info{1,15}='Pvֵ';dat_info{1,16}='Paֵ';dat_info{1,17}='Iaֵ';
 dat_info{1,18}='CAVֵ';dat_info{1,19}='IV2ֵ';dat_info{1,20}='aasֵ';dat_info{1,21}='avsֵ';
 dat_info{1,22}='adsֵ';
 %%%��ֵ�༰����������%%%%%
 
 %%%����������%%%%%%
 dat_info{1,23}='Tcֵ';dat_info{1,24}='TPֵ';dat_info{1,25}='Tvaֵ';
 %%%����������%%%%%%
 
 %%%����������%%%%%
 dat_info{1,26}='Ia_cֵ';dat_info{1,27}='Ia_bֵ';dat_info{1,28}='Bֵ';dat_info{1,29}='Cֵ';
 dat_info{1,30}='Pd_bֵ';dat_info{1,31}='Pd_cֵ';dat_info{1,32}='Pv_bֵ';dat_info{1,33}='Pv_cֵ';
 dat_info{1,34}='CAV_cֵ';dat_info{1,35}='CAV_bֵ';dat_info{1,36}='IV2_cֵ';dat_info{1,37}='IV2_bֵ';
 dat_info{1,38}='aas_cֵ';dat_info{1,39}='aas_bֵ';dat_info{1,40}='avs_cֵ';dat_info{1,41}='avs_bֵ';
 dat_info{1,42}='ads_cֵ';dat_info{1,43}='ads_bֵ';
 %%%����������%%%%%
 

 
 for f1 = 1:1:n1               %��1���ļ�����ѭ�������¼�
    ReN=strcat(cd_UD,w(f1).name,'\','*.UD');    %Ѱ�Ҽ�¼������
    ReC=dir(ReN);          %ͳ��һ�����ٸ���¼
    [n2,~]=size(ReC);      %�Լ�¼��������
    
    k=2;
    for f2 = 1:1:n2           %��2���ļ�����ѭ�������¼
        Rec_UD=strcat(cd_UD,w(f1).name,'\',ReC(f2).name);
        str_UD = Rec_UD;
        %----%%%%read data %%%%--
        %[DataUD,nt,t,ela,elo,sla,slo,M,sf,depth,Dir,Sta_H,Re,Maxacc,HD] = read_acc_data(str_UD);%Japan data read
        [DataUD,nt,t,ela,elo,sla,slo,M,Re,depth,site_con,M_type] = chinareaddata(str_UD);
        %Re���о�
        dt=1/nt;%ʱ������ntΪ����Ƶ��
        
        P_Acc = max(abs(DataUD));                  %��ֵ���ٶ�
        t_all=max(t);                              %��¼ʱ��
        
        Rec_EW=strcat(cd_EW,ReC(f2).name(1:end-3),".EW");
        str_EW=Rec_EW;
        %[DataEW,~,~] = read_acc_data_EW_NS(str_EW);
        [DataEW,~,~,~,~,~,~,~,~,~,~,~] = chinareaddata(str_EW);

        
        
        Rec_NS=strcat(cd_NS,ReC(f2).name(1:end-3),".NS");
        str_NS=Rec_NS;
        %[DataNS,nt,t] = read_acc_data_EW_NS(str_NS);
        [DataNS,nt,t,ela,elo,sla,slo,M,Re,depth,site_con,M_type] = chinareaddata(str_NS);
        
        
        %�˹��޸�P����ʱ
        [a1,a2] = textread('\P-wave arrival.txt','%s%f','headerlines',0);%a1�ļ�����a2 P����ʱ
        p_num=size(a1,1);
        for p_n_1=1:p_num
            if  strcmp(strcat(a1{p_n_1,1},'.UD'),ReC(f2).name)==1 %,'.UD'
                P_a_t_r=round(a2(p_n_1)*nt)/nt;%p����ʱ
                [SNR]= X_Z_B(DataUD,P_a_t_r,dt); %���������
                
                %%%��ֵ�༰����������%%%%%
                [Pd]= Pd_distance(DataUD,P_a_t_r,dt,3,4);%Pd
                [Pv]= Pv_distance(DataUD,P_a_t_r,dt,3,2);%Pv
                [Pa]= Pa_distance(DataUD,P_a_t_r,dt,3,2);%Pa
                [Ia]= Ia_distance(DataUD,DataEW,DataNS,P_a_t_r,dt,3,2);%Ia
                [CAV]= CAV_distance(DataUD,DataEW,DataNS,P_a_t_r,dt,3,2);%CAV
                [IV2]= IV2_distance(DataUD,DataEW,DataNS,P_a_t_r,dt,3,2);%IV2
                [aas]= aas_distance(DataUD,DataEW,DataNS,P_a_t_r,dt,3,4);%aas
                [avs]= avs_distance(DataUD,DataEW,DataNS,P_a_t_r,dt,3,4);%avs
                [ads]= ads_distance(DataUD,DataEW,DataNS,P_a_t_r,dt,3,4);%ads
                %%%��ֵ�༰����������%%%%%
                
                %%%����������%%%%%%
                [Tc,TP]= Tc_TP(DataUD,P_a_t_r,dt,3);%TC;TP
                [Tva]= Tva_feature(DataUD,P_a_t_r,dt,3);%Tva
                %disp;
                
                %%%����������%%%%%%
                
                %%%����������%%%%%
                [Ia_c]= Ia_c_distance(DataUD,DataEW,DataNS,P_a_t_r,dt,3,2); %������˹�Ҷ�ʱ��y=At
                [Ia_b]= Ia_b_distance(DataUD,DataEW,DataNS,P_a_t_r,dt,3,2);%������˹�Ҷ�ʱ��y=Bt*exp(-At)
                
                [B]= B_distance(DataUD,P_a_t_r,dt,3);%������ٶȰ���y=Bt*exp(-At) 
                [C]= C_distance(DataUD,P_a_t_r,dt,3,4);%������ٶȰ���y=Ct
                
                [Pd_b]= Pd_b_distance(DataUD,P_a_t_r,dt,3);%����λ�ư���y=Bt*exp(-At)
                [Pd_c]= Pd_c_distance(DataUD,P_a_t_r,dt,3,4);%����λ�ư���y=Ct
                
                [Pv_b]= Pv_b_distance(DataUD,P_a_t_r,dt,3);%�����ٶȰ���y=Bt*exp(-At)
                [Pv_c]= Pv_c_distance(DataUD,P_a_t_r,dt,3,4);%�����ٶȰ���y=Ct
                
                [CAV_c]= CAV_c_distance(DataUD,DataEW,DataNS,P_a_t_r,dt,3,2);%CAVʱ��y=C
                [CAV_b]= CAV_b_distance(DataUD,DataEW,DataNS,P_a_t_r,dt,3,2);%CAVʱ��y=Bt*exp(-At)
                
                [IV2_c]= IV2_c_distance(DataUD,DataEW,DataNS,P_a_t_r,dt,3,2);%IV2ʱ��y=Et
                [IV2_b]= IV2_b_distance(DataUD,DataEW,DataNS,P_a_t_r,dt,3,2);%IV2ʱ��y=Bt*exp(-At)
                
                [aas_c]= aas_c_distance(DataUD,DataEW,DataNS,P_a_t_r,dt,3,4);%aasʱ��y=Ct
                [aas_b]= aas_b_distance(DataUD,DataEW,DataNS,P_a_t_r,dt,3,4);%aasʱ��y=Bt*exp(-At)
                
                [avs_c]= avs_c_distance(DataUD,DataEW,DataNS,P_a_t_r,dt,3,4);%avsʱ��y=Ct
                [avs_b]= avs_b_distance(DataUD,DataEW,DataNS,P_a_t_r,dt,3,4);%avsʱ��y=Bt*exp(-At)
                
                [ads_c]= ads_c_distance(DataUD,DataEW,DataNS,P_a_t_r,dt,3,4);%adsʱ��y=Ct
                [ads_b]= ads_b_distance(DataUD,DataEW,DataNS,P_a_t_r,dt,3,4);%adsʱ��y=Bt*exp(-At)
                 %%%����������%%%%%
               
             
                

            end
        end
        

        
        
        dat_info{k,1}=SNR;dat_info{k,2}=P_a_t_r;dat_info{k,3}=M_type;
        dat_info{k,4}=ReC(f2).name;dat_info{k,5}=ela;dat_info{k,6}=elo;dat_info{k,7}=depth; 
        dat_info{k,8}=M;dat_info{k,9}=sla;dat_info{k,10}=slo;dat_info{k,11}=nt;
        dat_info{k,12}=Re;dat_info{k,13}=0;  
         
        %%%��ֵ�༰����������%%%%% 9
        dat_info{k,14}=log10(Pd);dat_info{k,15}=log10(Pv);dat_info{k,16}=log10(Pa);
        dat_info{k,17}=log10(Ia);dat_info{k,18}=log10(CAV);dat_info{k,19}=log10(IV2);
        dat_info{k,20}=log10(aas);dat_info{k,21}=log10(avs);dat_info{k,22}=log10(ads);
         %%%��ֵ�༰����������%%%%%
         
        %%%����������%%%%%% 3
        dat_info{k,23}=log10(Tc);dat_info{k,24}=log10(TP);dat_info{k,25}=log10(Tva);
        %%%����������%%%%%%
        
         %%%����������%%%%% 18
        dat_info{k,26}=log10(Ia_c);dat_info{k,27}=log10(Ia_b);dat_info{k,28}=log10(B);dat_info{k,29}=log10(C);
        dat_info{k,30}=log10(Pd_b);dat_info{k,31}=log10(Pd_c);dat_info{k,32}=log10(Pv_b);dat_info{k,33}=log10(Pv_c);
        dat_info{k,34}=log10(CAV_c);dat_info{k,35}=log10(CAV_b);dat_info{k,36}=log10(IV2_c);dat_info{k,37}=log10(IV2_b);
        dat_info{k,38}=log10(aas_c);dat_info{k,39}=log10(aas_b);dat_info{k,40}=log10(avs_c);dat_info{k,41}=log10(avs_b);
        dat_info{k,42}=log10(ads_c);dat_info{k,43}=log10(ads_b);
 %%%����������%%%%%
        
        movefile(Rec_UD,cd_move_UD);
        movefile(Rec_EW,cd_move_EW);
        movefile(Rec_NS,cd_move_NS);

        k=k+1;
    end
    
end
xlswrite('\GBR-Dis features.xlsx',dat_info);
