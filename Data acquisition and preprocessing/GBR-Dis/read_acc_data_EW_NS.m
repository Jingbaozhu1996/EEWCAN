function  [data,nt,t] = read_acc_data_EW_NS(str)
fid  = fopen(str);
%本程序目的：仅读取P波数据
%前提：1.顺序S/E/UP;2.每个数据集都是一样的；
%第1步，先读取整个数据的行数
%data-地震记录；ela、elo-地震纬度、经度；
%nt-采样频率；sla、slo-台站纬度、经度；
%t-持时；M-震级；sf-比例因子；Re-震中距；
%depth-震源深度；Sta_H-台站高度
%over = textscan(fid,'%s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s',inf)
%row  = length(over{:,1})
%第2步，读取头文件
%fclose(fid)
%fid  = fopen(str)
%over1 = textscan(fid,'%s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s',2*row/3)
line1 = textscan(fid,'%s %s %s %s %s %s %s %s %s %s %s %s %s %s %s',1);
line2 = textscan(fid,'%s %f ',1);
line3 = textscan(fid,'%s %f',1);
line4 = textscan(fid,'%s %s %f',1);
line5 = textscan(fid,'%s %f',1);
line6 = textscan(fid,'%s %s %s %s %s %s',1);
line7 = textscan(fid,'%s %s %f ',1);
line8 = textscan(fid,'%s %s %f',1);
line9 = textscan(fid,'%s %s %f ',1);
line10 = textscan(fid,'%s %s %s %s %s %s %s',1);
line11 = textscan(fid,'%s %s %f %s ',1);
line12 = textscan(fid,'%s %s %s',1);
line13 = textscan(fid,'%s %s',1);
line14 = textscan(fid,' %s %s %f(gal)/  %f ',1);
line15 = textscan(fid,' %s %s %s %s %s %s %s %s %s',1);
line16 = textscan(fid,' %s %s %s %s %s %s %s %s %s',1);
line17 = textscan(fid,' %s %s %s %s %s %s %s %s %s %s',1);
data    =  fscanf(fid,'%f',[1,inf]);
fclose(fid);
N  = length(data); %地震事件个数
nt = line11{3};          %采样频率
nt_s =1/nt;    %采样时间间隔
t = (0:1:N-1)*nt_s;            %持时
sla = line7{3};   %台站纬度
slo = line8{3};   %台站经度
ela = line2{2};   %地震纬度
elo = line3{2};   %地震经度
% Re=abs(111.12*cos(1/(sin(sla)*sin(ela)+cos(sla)*cos(ela)*cos(elo-slo))));%震中距
C = sin(sla/57.2958)*sin(ela/57.2958) + cos(sla/57.2958)*cos(ela/57.2958)*cos((slo-elo)/57.2958);
Re = 6371.004*acos(C);
sf = line14{3}/line14{4};
depth =line4{3};
M = line5{2};
Sta_H=line9{3};
Dir=char(line13{2});
Maxacc=line15{4};%最大加速度
sta_name = line6{3}; %台站名
HD=sqrt(depth^2+Re^2);%震源距
data=data*sf;
data = data - ones(size(data))*mean(data);   %基线矫正
UDf = max(abs(data));                                 %最大幅值