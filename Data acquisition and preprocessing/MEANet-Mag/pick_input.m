
clear;clc;
cd_UD ='\Data\';%数据在二级文件夹下
cd_EW ='\Data\EW\';%数据在一级文件夹下
cd_NS ='\Data\NS\';%数据在一级文件夹下

cd_move_UD = '\Data\pick\UD\';
cd_move_EW = '\Data\pick\EW\';
cd_move_NS = '\Data\pick\NS\';
mkdir(cd_move_UD);
mkdir(cd_move_EW);
mkdir(cd_move_NS);


g=dir(cd_UD);
w=g(3:end);                   %读取文件夹
[n1,~]=size(w);              %文件夹个数
 dat_info{1,1}='地震记录名称';dat_info{1,2}='震纬';dat_info{1,3}='震经';dat_info{1,4}='震源深度'; dat_info{1,5}='震级';
 dat_info{1,6}='台纬'; dat_info{1,7}='台经';dat_info{1,8}='采样频率';dat_info{1,9}='最大加速度';dat_info{1,10}='持时';
 dat_info{1,11}='震中距km';dat_info{1,12}='震源距km';dat_info{1,13}='信噪比';
 
 %xiang guan can shu
 PA_LVBO=zeros(10992,300);PD_LVBO=zeros(10992,300);PV_LVBO=zeros(10992,300);
 IV2_LVBO=zeros(10992,300);CAV_LVBO=zeros(10992,300);DE_LVBO=zeros(10992,300);
 cavv_LVBO=zeros(10992,300);caa_LVBO=zeros(10992,300);cad_LVBO=zeros(10992,300);
 tc_LVBO=zeros(10992,300);TP_LVBO=zeros(10992,300);Tva_LVBO=zeros(10992,300);
 ACC_LVBO=zeros(10992,300);VEL_LVBO=zeros(10992,300);DIS_LVBO=zeros(10992,300);
 OMG_ud=zeros(10992,300);
 
 
 z=1;
 k=2;
 for f1 = 1:1:n1               %在1级文件夹下循环地震事件
    ReN=strcat(cd_UD,w(f1).name,'\','*.UD');    %寻找记录的名字
    ReC=dir(ReN);           %统计一共多少个记录
    [n2,~]=size(ReC);      %对记录进行数数
    
    

    

    for f2 = 1:1:n2           %在2级文件夹下循环地震记录
        Rec_UD=strcat(cd_UD,w(f1).name,'\',ReC(f2).name);
        str_UD = Rec_UD;
        %----%%%%read data %%%%--
        [DataUD,nt,t,ela,elo,sla,slo,M,sf,depth,Dir,Sta_H,Re,Maxacc,HD] = read_acc_data(str_UD);
        %Re震中距
        dt=1/nt;%时间间隔，nt为采样频率
        
        P_Acc = max(abs(DataUD));                         %峰值加速度
        t_all=max(t);                                                %记录时长
        
        Rec_EW=strcat(cd_EW,ReC(f2).name(1:end-3),".EW");
        str_EW=Rec_EW;
        [DataEW,nt,t] = read_acc_data(str_EW);
        
        Rec_NS=strcat(cd_NS,ReC(f2).name(1:end-3),".NS");
        str_NS=Rec_NS;
        [DataNS,nt,t] = read_acc_data(str_NS);
        
        %人工修改P波到时
        [a1,a2] = textread('\P-wave arrival.txt','%s%f','headerlines',0);%a1文件名，a2 P波到时
        p_num=size(a1,1);
        for p_n_1=1:p_num
            if  strcmp(strcat(a1{p_n_1,1},'.UD'),ReC(f2).name)==1
                P_a_t_r=round(a2(p_n_1)*nt)/nt;%p波到时
 
                [SNR]= X_Z_B(DataUD,P_a_t_r,dt);
                [omg_ud]= omgud_lvbo_wave(DataUD,DataEW,DataNS,P_a_t_r,dt);%ud方向频率
                
                [pd_lvbo]= Pd_lvbo_wave(DataUD,DataEW,DataNS,P_a_t_r,dt); %%峰值位移波形
                [pa_lvbo,pv_lvbo]= Pa_Pv_lvbo_wave(DataUD,DataEW,DataNS,P_a_t_r,dt); %%
                [IV2_lvbo,Tva_lvbo,cavv_lvbo,caa_lvbo,DE_lvbo]= IV2_Tva_cav_caa_DE_lvbo_wave(DataUD,DataEW,DataNS,P_a_t_r,dt);%%
                [tc_lvbo,TP_lvbo,cad_lvbo]= tc_TP_cad_lvbo_wave(DataUD,DataEW,DataNS,P_a_t_r,dt);%%
                [CAV_lvbo]= CAV_lvbo_wave(DataUD,DataEW,DataNS,P_a_t_r,dt);%%
                [Dis_lvbo]= raw_lvbo_dis(DataUD,DataEW,DataNS,P_a_t_r,dt);%%
                [VEL_lvbo,ACC_lvbo]= raw_lvbo_acc_vel(DataUD,DataEW,DataNS,P_a_t_r,dt);%%
                %DISP;
                

                
            end
        end
        

        movefile(Rec_UD,cd_move_UD);
        movefile(Rec_EW,cd_move_EW);
        movefile(Rec_NS,cd_move_NS);
        
        
        PA_LVBO(z,:)=pa_lvbo;PD_LVBO(z,:)=pd_lvbo;PV_LVBO(z,:)=pv_lvbo;
        IV2_LVBO(z,:)=IV2_lvbo;CAV_LVBO(z,:)=CAV_lvbo; DE_LVBO(z,:)=DE_lvbo;
        cavv_LVBO(z,:)=cavv_lvbo;caa_LVBO(z,:)=caa_lvbo;cad_LVBO(z,:)=cad_lvbo;
        tc_LVBO(z,:)=tc_lvbo;TP_LVBO(z,:)=TP_lvbo; Tva_LVBO(z,:)=Tva_lvbo;
        ACC_LVBO(z,:)=ACC_lvbo;VEL_LVBO(z,:)=VEL_lvbo;DIS_LVBO(z,:)=Dis_lvbo;
        
        OMG_ud(z,:)=omg_ud;
        

        dat_info{k,1}=ReC(f2).name;dat_info{k,2}=ela;dat_info{k,3}=elo;dat_info{k,4}=depth; dat_info{k,5}=M;
        dat_info{k,6}=sla;dat_info{k,7}=slo;dat_info{k,8}=nt;dat_info{k,9}=P_Acc;dat_info{k,10}=t_all;
        dat_info{k,11}=Re;dat_info{k,12}=HD; dat_info{k,13}=SNR; 

       
        k=k+1;
        z=z+1;

        % pause
    end
    
end
xlswrite('\MEANet-Mag\xinxi_3s.xlsx',dat_info);


save('\MEANet-Mag\PA_LVBO.mat','PA_LVBO');
save('\MEANet-Mag\PD_LVBO.mat','PD_LVBO');
save('\MEANet-Mag\PV_LVBO.mat','PV_LVBO');
save('\MEANet-Mag\IV2_LVBO.mat','IV2_LVBO');
save('\MEANet-Mag\CAV_LVBO.mat','CAV_LVBO');
save('\MEANet-Mag\DE_LVBO.mat','DE_LVBO');
save('\MEANet-Mag\cavv_LVBO.mat','cavv_LVBO');
save('\MEANet-Mag\caa_LVBO.mat','caa_LVBO');
save('\MEANet-Mag\cad_LVBO.mat','cad_LVBO');
save('\MEANet-Mag\tc_LVBO.mat','tc_LVBO');
save('\MEANet-Mag\TP_LVBO.mat','TP_LVBO');
save('\MEANet-Mag\Tva_LVBO.mat','Tva_LVBO');
save('\MEANet-Mag\ACC_LVBO.mat','ACC_LVBO');
save('\MEANet-Mag\VEL_LVBO.mat','VEL_LVBO');
save('\MEANet-Mag\DIS_LVBO.mat','DIS_LVBO');
save('\MEANet-Mag\OMG_ud.mat','OMG_ud');