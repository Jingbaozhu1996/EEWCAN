function  [data,nt,t] = read_acc_data_EW_NS(str)
fid  = fopen(str);
%������Ŀ�ģ�����ȡP������
%ǰ�᣺1.˳��S/E/UP;2.ÿ�����ݼ�����һ���ģ�
%��1�����ȶ�ȡ�������ݵ�����
%data-�����¼��ela��elo-����γ�ȡ����ȣ�
%nt-����Ƶ�ʣ�sla��slo-̨վγ�ȡ����ȣ�
%t-��ʱ��M-�𼶣�sf-�������ӣ�Re-���оࣻ
%depth-��Դ��ȣ�Sta_H-̨վ�߶�
%over = textscan(fid,'%s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s',inf)
%row  = length(over{:,1})
%��2������ȡͷ�ļ�
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
N  = length(data); %�����¼�����
nt = line11{3};          %����Ƶ��
nt_s =1/nt;    %����ʱ����
t = (0:1:N-1)*nt_s;            %��ʱ
sla = line7{3};   %̨վγ��
slo = line8{3};   %̨վ����
ela = line2{2};   %����γ��
elo = line3{2};   %���𾭶�
% Re=abs(111.12*cos(1/(sin(sla)*sin(ela)+cos(sla)*cos(ela)*cos(elo-slo))));%���о�
C = sin(sla/57.2958)*sin(ela/57.2958) + cos(sla/57.2958)*cos(ela/57.2958)*cos((slo-elo)/57.2958);
Re = 6371.004*acos(C);
sf = line14{3}/line14{4};
depth =line4{3};
M = line5{2};
Sta_H=line9{3};
Dir=char(line13{2});
Maxacc=line15{4};%�����ٶ�
sta_name = line6{3}; %̨վ��
HD=sqrt(depth^2+Re^2);%��Դ��
data=data*sf;
data = data - ones(size(data))*mean(data);   %���߽���
UDf = max(abs(data));                                 %����ֵ