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
w=g(3:end);                   %读取文件夹
[n1,~]=size(w);              %文件夹个数

 %xiang guan can shu
 dat_info{1,1}='峰值位移Pd';dat_info{1,2}='峰值速度Pv';dat_info{1,3}='峰值加速度Pa';
 dat_info{1,4}='速度平方积分IV2';dat_info{1,5}='cav单分向ud的累积';dat_info{1,6}='caa单分向ud的累积';
 dat_info{1,7}='cad单分向ud的累积';dat_info{1,8}='CAV累积绝对速度(三分向)';
 dat_info{1,9}='累积能量变化率DE';dat_info{1,10}='构造参数tc*Pd';dat_info{1,11}='特征周期tc';
 dat_info{1,12}='峰值比Tva'; dat_info{1,13}='震级';dat_info{1,14}='震中距km';dat_info{1,15}='震源距km';
 dat_info{1,16}='震纬';dat_info{1,17}='震经';dat_info{1,18}='震源深度';dat_info{1,19}='信噪比';
 
 dat_info{1,20}='地震记录名称';dat_info{1,21}='台纬'; dat_info{1,22}='台经';dat_info{1,23}='采样频率';
 dat_info{1,24}='最大加速度';dat_info{1,25}='持时';dat_info{1,26}='P波到时'; 
 

 for f1 = 1:1:n1               %在1级文件夹下循环地震事件
    ReN=strcat(cd_UD,w(f1).name,'\','*.UD');    %寻找记录的名字
    ReC=dir(ReN);           %统计一共多少个记录
    [n2,~]=size(ReC);      %对记录进行数数
    
    k=2;
    for f2 = 1:1:n2           %在2级文件夹下循环地震记录
        Rec_UD=strcat(cd_UD,w(f1).name,'\',ReC(f2).name);
        str_UD = Rec_UD;
        %----%%%%read data %%%%--
        [DataUD,nt,t,ela,elo,sla,slo,M,sf,depth,Dir,Sta_H,Re,Maxacc,HD] = read_acc_data(str_UD);
        %Re震中距
        dt=1/nt; %时间间隔，nt为采样频率
        
        P_Acc = max(abs(DataUD));                         %峰值加速度
        t_all=max(t);                                                %记录时长
        
        Rec_EW=strcat(cd_EW,ReC(f2).name(1:end-3),".EW");
        str_EW=Rec_EW;
        [DataEW,nt,t] = read_acc_data_EW_NS(str_EW);
        
        Rec_NS=strcat(cd_NS,ReC(f2).name(1:end-3),".NS");
        str_NS=Rec_NS;
        [DataNS,nt,t] = read_acc_data_EW_NS(str_NS);
        
        %人工修改P波到时
        [a1,a2] = textread('\P-wave arrival.txt','%s%f','headerlines',0);%a1文件名，a2 P波到时
        p_num=size(a1,1);
        for p_n_1=1:p_num
            if  strcmp(strcat(a1{p_n_1,1},'.UD'),ReC(f2).name)==1 %,'.UD'
                P_a_t_r=round(a2(p_n_1)*nt)/nt;%p波到时 
                [SNR] = X_Z_B(DataUD,P_a_t_r,dt);
                [tc,Pd,TP,cad]= xzb_d_cs(DataUD,P_a_t_r,dt);
           
                [Pv,Pa,IV2,Tva,cav,caa,DE]= xzb_x_cs(DataUD,P_a_t_r,dt);
                [a_ud,t2]= bate_lvbo(DataUD,P_a_t_r,dt);%滤波后，时间窗内的加速度,t2是时间窗,UDfangxiang
                [a_ew,t2]= bate_lvbo(DataEW,P_a_t_r,dt);%滤波后，时间窗内的加速度,EWfangxiang
                [a_ns,t2]= bate_lvbo(DataNS,P_a_t_r,dt);%滤波后，时间窗内的加速度,NSfangxiang
                CAA= sqrt(((a_ud).^2)+((a_ew).^2)+((a_ns).^2));%时间窗内三分向合成加速度
                CAV=max(cumtrapz(t2,abs(CAA)));%累积绝对速度，三分向的
                
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
