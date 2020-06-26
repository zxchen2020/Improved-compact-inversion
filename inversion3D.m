clear all;
x0=input('Enter the Model starting coordinates x0');
y0=input('Enter the Model starting coordinates y0');
z0=input('Enter the Model starting coordinates z0');%x0,y0,z0Ϊģ����ʼλ�õ�����
Mx=input('Enter the Number of model grids Mx:');
My=input('Enter the Number of model grids My:');
Mz=input('Enter the Number of model grids Mz:');%Mx,My,MzΪģ����������xΪ����yΪ����
Dx=input('Enter the Model grid size Dx:');
Dy=input('Enter the Model grid size Dy:');
Dz=input('Enter the Model grid size Dz:');%Dx,Dy,DzΪģ�������С(m)
l_0=input('Enter the SNR l_0:');%ʵ�����������
bet=input('Enter the weight coefficient bet:');%ģ��Լ�������Ȩ��ϵ��
Max_iter=input('Enter the Maximum number of iterations Max_iter:');%����������
Max_p=input ('Enter the Upper density limit Max_p');%�ܶ�����ֵ
Min_p=input ('Enter the Upper density limit Min_p');%�ܶ�����ֵ
%����Ĳ���
M=Mx*My*Mz;
d_obs=load('d_obs.txt');%ʵ�������������ݶ�����
N=size(d_obs,1);
P=d_obs(:,1:3);%ʵ���쳣�۲������
dxy=d_obs(:,6)*10^(-9);%���ݼ����������ݣ�d_obs�ĵ�λ�ֱ�ΪmGal��E�������ݼ���ʱ���õ�λ�ֱ�Ϊm/s^2=10^5mGal�� s^-2=10^9E
[Gz,Gxx,Gxy,Gxz,Gyy,Gyz,Gzz]=Grav3D(N,P,M,x0,y0,z0,Mx,My,Mz,Dx,Dy,Dz);%����3D���˾���
A=Gxy;d=dxy;
W_xyz=Matrix_smooth(M,Mx,My,Mz);%����⻬Լ������
I_n=eye(N);         
I_m=eye(M);
W_z=zeros(M,M);%��ȼ�Ȩ����
p_old=zeros(M,1);p_new=zeros(M,1);
c=sum(A.^(2));
for i=1:M   
    W_z(i,i)=c(i)^(-0.25);
end
Wm_0=I_m;%��ʼ����Լ������
Wm=(1-bet)*Wm_0+bet*W_xyz;%��ʼ������ģ��Լ������
D=(A/Wm*W_z*W_z*A').*I_n;
W_e=l_0^(2)*D;%��ʼ����������Ȩ����   
eps_sus=[10^(-10),10^(-11),10^(-12),10^(-13)];%�Ŷ�ϵ��                              
for l=1:length(eps_sus)              
    for k=1:Max_iter 
        p_old=p_new;
        p_new=Wm\W_z*W_z*A'/(A/Wm*W_z*W_z*A'+W_e)*d;%���ݼ��㣬p_newΪ��������ģ�͵��ܶȷֲ�(g/cm^3)
        Wm_0=(((diag(p_new)).^(2))+(eps_sus(l))*I_m)^-1;%�������Լ������
        Wm=(1-bet)*Wm_0+bet*W_xyz;%������Լ������
        D=(A/Wm*W_z*W_z*A').*I_n;
        W_e=l_0^(2)*D;%����������Ȩ����
        min=find(p_new<Min_p);  
        max=find(p_new>Max_p);
        p_new(min)=Min_p;                
        p_new(max)=Max_p;%����Լ��
        if abs(p_new-p_old)<(0.01*abs(p_old))
            break;%ֹͣ�����ı�׼��ģ�͸ı���
        end    
     end           
end
    
    
      