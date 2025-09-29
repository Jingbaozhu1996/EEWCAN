function    [ data4,nt,t,ela,elo,sla,slo,M,R,depth,site_con,M_type] =chinareaddata(str2)
        %UD
        fid1=fopen(str2);
        line1=textscan(fid1,'%s',1);
        line2=textscan(fid1,'%s %s %s %s',1);                       % 发震时间
        line3=textscan(fid1, '%[^\n]',1);                                % 地震地点
        line4_1=textscan(fid1,'EPICENTER  %fN',1);                              
        ela=line4_1{1};                                                           % 震中纬度
        line4_2=textscan(fid1,'%fE',1);                                    % 震中经度
        elo=line4_2{1};
        %line4_3=textscan(fid1, 'DEPTH %f KM ',1);                   % 震源深度
        line4_3=textscan(fid1, '%s %s %s ',1);       % 震源深度
        line5=textscan(fid1,'MAGNITUDE %f(%s) ',1);              % 震级
        line6=textscan(fid1,'%s %s %fN %fE',1);                      % 台站编号
        sla=line6{3};
        slo=line6{4};
        line7=textscan(fid1,'SITE CONDITION: %s',1);              % 场地条件
        line8=textscan(fid1,'INSTRUMENT TYPE: %s',1);           % 仪器型号
        line9=textscan(fid1,'OBSERVING POINT: %s',1);           % 观测点位置
        line10=textscan(fid1,'%s %s',1);                                   % 方向
        line11_1=textscan(fid1,'%s %s',1);                                % 未校正记录
        line11_2=textscan(fid1,'UNIT: %s',1);                            % 单位
        line12_1=textscan(fid1,'NO. OF POINTS:  %d',1);           % 采样点数
        line12_2=textscan(fid1,'EQUALLY SPACED INTERVALS OF:  %f  SEC',1);      % 时间间隔
        line13_1=textscan(fid1,'PEAK VALUE:    %f   AT',1);        % 峰值
        line13_2=textscan(fid1,'%f SEC',1);                                % 峰值时间
        line13_3=textscan(fid1,'DURATION: %d SEC',1);             % 持时
        line14=textscan(fid1, '%[^\n]',1); 
        line15=textscan(fid1,'%s',1);                                           % 读取CSMNC
        data2=fscanf(fid1,'%f',[1,inf]);                                        %读数据时是按行逐个读取,排成一行(读取数据，读成一行)
        fclose(fid1);
        %地震记录不从零开始，去均值减掉全部记录的平均,（基线调整）
        N=length(data2);                                                         %求地震记录的个数
        nt=double(1/double(line12_2{1}));                                %采样频率（化成数）（元胞数组化成数）
        t=(1:N)*(double(line12_2{1}));                                        %时间间隔是line12_2{1},x轴是时间t
        M=line5{1};                                                                   %震级
        M_type=line5{2}{1,1}(1:2);                                              %震级
        AD2=ones(size(data2))*mean(data2);                             %求全部地震记录的平均
        data4=data2-AD2;                                                         %基线校正_真实记录
        UDf=max(abs(data4));    
        %计算震中距
        Dis=distance(ela,elo,sla,slo);
        R=Dis*6371*2*pi/360;
        depth = line4_3{1,2}{1,1};
        site_con=line7{1,1}{1,1};
end