function    [ data4,nt,t,ela,elo,sla,slo,M,R,depth,site_con,M_type] =chinareaddata(str2)
        %UD
        fid1=fopen(str2);
        line1=textscan(fid1,'%s',1);
        line2=textscan(fid1,'%s %s %s %s',1);                       % ����ʱ��
        line3=textscan(fid1, '%[^\n]',1);                                % ����ص�
        line4_1=textscan(fid1,'EPICENTER  %fN',1);                              
        ela=line4_1{1};                                                           % ����γ��
        line4_2=textscan(fid1,'%fE',1);                                    % ���о���
        elo=line4_2{1};
        %line4_3=textscan(fid1, 'DEPTH %f KM ',1);                   % ��Դ���
        line4_3=textscan(fid1, '%s %s %s ',1);       % ��Դ���
        line5=textscan(fid1,'MAGNITUDE %f(%s) ',1);              % ��
        line6=textscan(fid1,'%s %s %fN %fE',1);                      % ̨վ���
        sla=line6{3};
        slo=line6{4};
        line7=textscan(fid1,'SITE CONDITION: %s',1);              % ��������
        line8=textscan(fid1,'INSTRUMENT TYPE: %s',1);           % �����ͺ�
        line9=textscan(fid1,'OBSERVING POINT: %s',1);           % �۲��λ��
        line10=textscan(fid1,'%s %s',1);                                   % ����
        line11_1=textscan(fid1,'%s %s',1);                                % δУ����¼
        line11_2=textscan(fid1,'UNIT: %s',1);                            % ��λ
        line12_1=textscan(fid1,'NO. OF POINTS:  %d',1);           % ��������
        line12_2=textscan(fid1,'EQUALLY SPACED INTERVALS OF:  %f  SEC',1);      % ʱ����
        line13_1=textscan(fid1,'PEAK VALUE:    %f   AT',1);        % ��ֵ
        line13_2=textscan(fid1,'%f SEC',1);                                % ��ֵʱ��
        line13_3=textscan(fid1,'DURATION: %d SEC',1);             % ��ʱ
        line14=textscan(fid1, '%[^\n]',1); 
        line15=textscan(fid1,'%s',1);                                           % ��ȡCSMNC
        data2=fscanf(fid1,'%f',[1,inf]);                                        %������ʱ�ǰ��������ȡ,�ų�һ��(��ȡ���ݣ�����һ��)
        fclose(fid1);
        %�����¼�����㿪ʼ��ȥ��ֵ����ȫ����¼��ƽ��,�����ߵ�����
        N=length(data2);                                                         %������¼�ĸ���
        nt=double(1/double(line12_2{1}));                                %����Ƶ�ʣ�����������Ԫ�����黯������
        t=(1:N)*(double(line12_2{1}));                                        %ʱ������line12_2{1},x����ʱ��t
        M=line5{1};                                                                   %��
        M_type=line5{2}{1,1}(1:2);                                              %��
        AD2=ones(size(data2))*mean(data2);                             %��ȫ�������¼��ƽ��
        data4=data2-AD2;                                                         %����У��_��ʵ��¼
        UDf=max(abs(data4));    
        %�������о�
        Dis=distance(ela,elo,sla,slo);
        R=Dis*6371*2*pi/360;
        depth = line4_3{1,2}{1,1};
        site_con=line7{1,1}{1,1};
end