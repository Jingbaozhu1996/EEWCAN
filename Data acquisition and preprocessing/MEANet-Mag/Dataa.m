
function [a1,d1,v1]= Dataa(inputA,output2,dt)

[m n]=size(inputA);%输入加速度的矩阵维度
for j=1:n
    t(j,1)=(j-1)*dt;
end
t4=3;%表示时间窗长度s
v=cumtrapz(t,inputA);%速度，这是对加速度求积分
d=cumtrapz(t,v);%位移，这是对速度求积分
[~,q]=min(abs(t-output2));
% t1=output2:output2+3;
n1=find(t==t(q,1));%n1是P波到时的后一个点
v1=abs(v(n1+1:n1+t4/dt));%p波到时后时间窗内的速度幅值
d1=abs(d(n1+1:n1+t4/dt));%时间窗内的位移
a1=abs(inputA(n1+1:n1+t4/dt));%P波后3s的加速度
if n1-2-t4/dt<=0
    v2=max(abs(v(1:n1-2)));%P波前3s的速度幅值
    a2=inputA(1:n1-2);
else
    v2=max(abs(v(n1-2-t4/dt:n1-2)));%P波前3s的速度幅值
    a2=inputA(n1-2-t4/dt:n1-2);
end
SNR=v1/v2;%信噪比
end

