clear all;
x0=input('Enter the Model starting coordinates x0');
y0=input('Enter the Model starting coordinates y0');
z0=input('Enter the Model starting coordinates z0');%x0,y0,z0为模型起始位置的坐标
Mx=input('Enter the Number of model grids Mx:');
My=input('Enter the Number of model grids My:');
Mz=input('Enter the Number of model grids Mz:');%Mx,My,Mz为模型网格数，x为西向，y为北向
Dx=input('Enter the Model grid size Dx:');
Dy=input('Enter the Model grid size Dy:');
Dz=input('Enter the Model grid size Dz:');%Dx,Dy,Dz为模型网格大小(m)
l_0=input('Enter the SNR l_0:');%实测数据信噪比
bet=input('Enter the weight coefficient bet:');%模型约束矩阵的权重系数
Max_iter=input('Enter the Maximum number of iterations Max_iter:');%最大迭代次数
Max_p=input ('Enter the Upper density limit Max_p');%密度上限值
Min_p=input ('Enter the Upper density limit Min_p');%密度下限值
%输入的参数
M=Mx*My*Mz;
d_obs=load('d_obs.txt');%实测重力及重力梯度数据
N=size(d_obs,1);
P=d_obs(:,1:3);%实测异常观测点坐标
dxy=d_obs(:,6)*10^(-9);%反演计算所以数据，d_obs的单位分别为mGal和E，而反演计算时所用单位分别为m/s^2=10^5mGal、 s^-2=10^9E
[Gz,Gxx,Gxy,Gxz,Gyy,Gyz,Gzz]=Grav3D(N,P,M,x0,y0,z0,Mx,My,Mz,Dx,Dy,Dz);%计算3D主核矩阵
A=Gxy;d=dxy;
W_xyz=Matrix_smooth(M,Mx,My,Mz);%计算光滑约束矩阵
I_n=eye(N);         
I_m=eye(M);
W_z=zeros(M,M);%深度加权矩阵
p_old=zeros(M,1);p_new=zeros(M,1);
c=sum(A.^(2));
for i=1:M   
    W_z(i,i)=c(i)^(-0.25);
end
Wm_0=I_m;%初始紧凑约束矩阵
Wm=(1-bet)*Wm_0+bet*W_xyz;%初始计算混合模型约束矩阵
D=(A/Wm*W_z*W_z*A').*I_n;
W_e=l_0^(2)*D;%初始计算噪声加权矩阵   
eps_sus=[10^(-10),10^(-11),10^(-12),10^(-13)];%扰动系数                              
for l=1:length(eps_sus)              
    for k=1:Max_iter 
        p_old=p_new;
        p_new=Wm\W_z*W_z*A'/(A/Wm*W_z*W_z*A'+W_e)*d;%反演计算，p_new为计算所得模型的密度分布(g/cm^3)
        Wm_0=(((diag(p_new)).^(2))+(eps_sus(l))*I_m)^-1;%计算紧凑约束矩阵
        Wm=(1-bet)*Wm_0+bet*W_xyz;%计算混合约束矩阵
        D=(A/Wm*W_z*W_z*A').*I_n;
        W_e=l_0^(2)*D;%计算噪声加权矩阵
        min=find(p_new<Min_p);  
        max=find(p_new>Max_p);
        p_new(min)=Min_p;                
        p_new(max)=Max_p;%物性约束
        if abs(p_new-p_old)<(0.01*abs(p_old))
            break;%停止迭代的标准，模型改变量
        end    
     end           
end
    
    
      