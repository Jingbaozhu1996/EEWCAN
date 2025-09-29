clear;clc;
cd_UD = '\Data\';%数据在二级文件夹下
cd_EW = '\Data\EW\';%数据在一级文件夹下
cd_NS = '\Data\NS\';%数据在一级文件夹下

cd_move_UD = '\Data\pick\UD\';
cd_move_EW = '\Data\pick\EW\';
cd_move_NS = '\Data\pick\NS\';

mkdir(cd_move_UD);
mkdir(cd_move_EW);
mkdir(cd_move_NS);



g=dir(cd_UD);
w=g(3:end);                  %读取文件夹
[n1,~]=size(w);              %文件夹个数

 dat_info{1,1}='信噪比';dat_info{1,2}='P波到时';dat_info{1,3}='M_type';
 dat_info{1,4}='地震记录名称';dat_info{1,5}='震纬';dat_info{1,6}='震经';dat_info{1,7}='震源深度'; 
 dat_info{1,8}='震级';dat_info{1,9}='台纬'; dat_info{1,10}='台经';dat_info{1,11}='采样频率';
 dat_info{1,12}='震中距km';dat_info{1,13}='震源距km';
 
 %%%幅值类及能量类特征%%%%%
 dat_info{1,14}='Pd值';dat_info{1,15}='Pv值';dat_info{1,16}='Pa值';dat_info{1,17}='Ia值';
 dat_info{1,18}='CAV值';dat_info{1,19}='IV2值';dat_info{1,20}='aas值';dat_info{1,21}='avs值';
 dat_info{1,22}='ads值';
 %%%幅值类及能量类特征%%%%%
 
 %%%周期类特征%%%%%%
 dat_info{1,23}='Tc值';dat_info{1,24}='TP值';dat_info{1,25}='Tva值';
 %%%周期类特征%%%%%%
 
 %%%包络增长率%%%%%
 dat_info{1,26}='Ia_c值';dat_info{1,27}='Ia_b值';dat_info{1,28}='B值';dat_info{1,29}='C值';
 dat_info{1,30}='Pd_b值';dat_info{1,31}='Pd_c值';dat_info{1,32}='Pv_b值';dat_info{1,33}='Pv_c值';
 dat_info{1,34}='CAV_c值';dat_info{1,35}='CAV_b值';dat_info{1,36}='IV2_c值';dat_info{1,37}='IV2_b值';
 dat_info{1,38}='aas_c值';dat_info{1,39}='aas_b值';dat_info{1,40}='avs_c值';dat_info{1,41}='avs_b值';
 dat_info{1,42}='ads_c值';dat_info{1,43}='ads_b值';
 %%%包络增长率%%%%%
 

 
 for f1 = 1:1:n1               %在1级文件夹下循环地震事件
    ReN=strcat(cd_UD,w(f1).name,'\','*.UD');    %寻找记录的名字
    ReC=dir(ReN);          %统计一共多少个记录
    [n2,~]=size(ReC);      %对记录进行数数
    
    k=2;
    for f2 = 1:1:n2           %在2级文件夹下循环地震记录
        Rec_UD=strcat(cd_UD,w(f1).name,'\',ReC(f2).name);
        str_UD = Rec_UD;
        %----%%%%read data %%%%--
        %[DataUD,nt,t,ela,elo,sla,slo,M,sf,depth,Dir,Sta_H,Re,Maxacc,HD] = read_acc_data(str_UD);%Japan data read
        [DataUD,nt,t,ela,elo,sla,slo,M,Re,depth,site_con,M_type] = chinareaddata(str_UD);
        %Re震中距
        dt=1/nt;%时间间隔，nt为采样频率
        
        P_Acc = max(abs(DataUD));                  %峰值加速度
        t_all=max(t);                              %记录时长
        
        Rec_EW=strcat(cd_EW,ReC(f2).name(1:end-3),".EW");
        str_EW=Rec_EW;
        %[DataEW,~,~] = read_acc_data_EW_NS(str_EW);
        [DataEW,~,~,~,~,~,~,~,~,~,~,~] = chinareaddata(str_EW);

        
        
        Rec_NS=strcat(cd_NS,ReC(f2).name(1:end-3),".NS");
        str_NS=Rec_NS;
        %[DataNS,nt,t] = read_acc_data_EW_NS(str_NS);
        [DataNS,nt,t,ela,elo,sla,slo,M,Re,depth,site_con,M_type] = chinareaddata(str_NS);
        
        
        %人工修改P波到时
        [a1,a2] = textread('\P-wave arrival.txt','%s%f','headerlines',0);%a1文件名，a2 P波到时
        p_num=size(a1,1);
        for p_n_1=1:p_num
            if  strcmp(strcat(a1{p_n_1,1},'.UD'),ReC(f2).name)==1 %,'.UD'
                P_a_t_r=round(a2(p_n_1)*nt)/nt;%p波到时
                [SNR]= X_Z_B(DataUD,P_a_t_r,dt); %计算信噪比
                
                %%%幅值类及能量类特征%%%%%
                [Pd]= Pd_distance(DataUD,P_a_t_r,dt,3,4);%Pd
                [Pv]= Pv_distance(DataUD,P_a_t_r,dt,3,2);%Pv
                [Pa]= Pa_distance(DataUD,P_a_t_r,dt,3,2);%Pa
                [Ia]= Ia_distance(DataUD,DataEW,DataNS,P_a_t_r,dt,3,2);%Ia
                [CAV]= CAV_distance(DataUD,DataEW,DataNS,P_a_t_r,dt,3,2);%CAV
                [IV2]= IV2_distance(DataUD,DataEW,DataNS,P_a_t_r,dt,3,2);%IV2
                [aas]= aas_distance(DataUD,DataEW,DataNS,P_a_t_r,dt,3,4);%aas
                [avs]= avs_distance(DataUD,DataEW,DataNS,P_a_t_r,dt,3,4);%avs
                [ads]= ads_distance(DataUD,DataEW,DataNS,P_a_t_r,dt,3,4);%ads
                %%%幅值类及能量类特征%%%%%
                
                %%%周期类特征%%%%%%
                [Tc,TP]= Tc_TP(DataUD,P_a_t_r,dt,3);%TC;TP
                [Tva]= Tva_feature(DataUD,P_a_t_r,dt,3);%Tva
                %disp;
                
                %%%周期类特征%%%%%%
                
                %%%包络增长率%%%%%
                [Ia_c]= Ia_c_distance(DataUD,DataEW,DataNS,P_a_t_r,dt,3,2); %阿里亚斯烈度时程y=At
                [Ia_b]= Ia_b_distance(DataUD,DataEW,DataNS,P_a_t_r,dt,3,2);%阿里亚斯烈度时程y=Bt*exp(-At)
                
                [B]= B_distance(DataUD,P_a_t_r,dt,3);%竖向加速度包络y=Bt*exp(-At) 
                [C]= C_distance(DataUD,P_a_t_r,dt,3,4);%竖向加速度包络y=Ct
                
                [Pd_b]= Pd_b_distance(DataUD,P_a_t_r,dt,3);%竖向位移包络y=Bt*exp(-At)
                [Pd_c]= Pd_c_distance(DataUD,P_a_t_r,dt,3,4);%竖向位移包络y=Ct
                
                [Pv_b]= Pv_b_distance(DataUD,P_a_t_r,dt,3);%竖向速度包络y=Bt*exp(-At)
                [Pv_c]= Pv_c_distance(DataUD,P_a_t_r,dt,3,4);%竖向速度包络y=Ct
                
                [CAV_c]= CAV_c_distance(DataUD,DataEW,DataNS,P_a_t_r,dt,3,2);%CAV时程y=C
                [CAV_b]= CAV_b_distance(DataUD,DataEW,DataNS,P_a_t_r,dt,3,2);%CAV时程y=Bt*exp(-At)
                
                [IV2_c]= IV2_c_distance(DataUD,DataEW,DataNS,P_a_t_r,dt,3,2);%IV2时程y=Et
                [IV2_b]= IV2_b_distance(DataUD,DataEW,DataNS,P_a_t_r,dt,3,2);%IV2时程y=Bt*exp(-At)
                
                [aas_c]= aas_c_distance(DataUD,DataEW,DataNS,P_a_t_r,dt,3,4);%aas时程y=Ct
                [aas_b]= aas_b_distance(DataUD,DataEW,DataNS,P_a_t_r,dt,3,4);%aas时程y=Bt*exp(-At)
                
                [avs_c]= avs_c_distance(DataUD,DataEW,DataNS,P_a_t_r,dt,3,4);%avs时程y=Ct
                [avs_b]= avs_b_distance(DataUD,DataEW,DataNS,P_a_t_r,dt,3,4);%avs时程y=Bt*exp(-At)
                
                [ads_c]= ads_c_distance(DataUD,DataEW,DataNS,P_a_t_r,dt,3,4);%ads时程y=Ct
                [ads_b]= ads_b_distance(DataUD,DataEW,DataNS,P_a_t_r,dt,3,4);%ads时程y=Bt*exp(-At)
                 %%%包络增长率%%%%%
               
             
                

            end
        end
        

        
        
        dat_info{k,1}=SNR;dat_info{k,2}=P_a_t_r;dat_info{k,3}=M_type;
        dat_info{k,4}=ReC(f2).name;dat_info{k,5}=ela;dat_info{k,6}=elo;dat_info{k,7}=depth; 
        dat_info{k,8}=M;dat_info{k,9}=sla;dat_info{k,10}=slo;dat_info{k,11}=nt;
        dat_info{k,12}=Re;dat_info{k,13}=0;  
         
        %%%幅值类及能量类特征%%%%% 9
        dat_info{k,14}=log10(Pd);dat_info{k,15}=log10(Pv);dat_info{k,16}=log10(Pa);
        dat_info{k,17}=log10(Ia);dat_info{k,18}=log10(CAV);dat_info{k,19}=log10(IV2);
        dat_info{k,20}=log10(aas);dat_info{k,21}=log10(avs);dat_info{k,22}=log10(ads);
         %%%幅值类及能量类特征%%%%%
         
        %%%周期类特征%%%%%% 3
        dat_info{k,23}=log10(Tc);dat_info{k,24}=log10(TP);dat_info{k,25}=log10(Tva);
        %%%周期类特征%%%%%%
        
         %%%包络增长率%%%%% 18
        dat_info{k,26}=log10(Ia_c);dat_info{k,27}=log10(Ia_b);dat_info{k,28}=log10(B);dat_info{k,29}=log10(C);
        dat_info{k,30}=log10(Pd_b);dat_info{k,31}=log10(Pd_c);dat_info{k,32}=log10(Pv_b);dat_info{k,33}=log10(Pv_c);
        dat_info{k,34}=log10(CAV_c);dat_info{k,35}=log10(CAV_b);dat_info{k,36}=log10(IV2_c);dat_info{k,37}=log10(IV2_b);
        dat_info{k,38}=log10(aas_c);dat_info{k,39}=log10(aas_b);dat_info{k,40}=log10(avs_c);dat_info{k,41}=log10(avs_b);
        dat_info{k,42}=log10(ads_c);dat_info{k,43}=log10(ads_b);
 %%%包络增长率%%%%%
        
        movefile(Rec_UD,cd_move_UD);
        movefile(Rec_EW,cd_move_EW);
        movefile(Rec_NS,cd_move_NS);

        k=k+1;
    end
    
end
xlswrite('\GBR-Dis features.xlsx',dat_info);
