
function [Gz,Gxx,Gxy,Gxz,Gyy,Gyz,Gzz]=Grav3D(N,P,M,x0,y0,z0,Mx,My,Mz,Dx,Dy,Dz)
K=zeros(M,3);%K是各个长方体块模型的中心点坐标
Gz=zeros(N,M);
Gxx=zeros(N,M);
Gxy=zeros(N,M);
Gxz=zeros(N,M);
Gyy=zeros(N,M);
Gyz=zeros(N,M);
Gzz=zeros(N,M);%G主核矩阵
for i=1:My
    yk=y0+Dy/2+(i-1)*Dy;
    for j=1:Mx
        xk=x0+Dx/2+(j-1)*Dx;
        for k=1:Mz
        zk=z0+Dz/2+(k-1)*Dz;
        l=k+(i-1)*Mx*Mz+(j-1)*Mz;
        K(l,1)=xk;
        K(l,2)=yk;
        K(l,3)=zk;
        end
    end
end
for n=1:N
    for m=1:M
        a=6.67259*10^-8;%引力常量cm^3/(g.s^2)%
        Gz(n,m)=a*((K(m,3)-P(n,3))*Dx*Dy*Dz)/((K(m,1)-P(n,1))^2+(K(m,2)-P(n,2))^2+(K(m,3)-P(n,3))^2)^(3/2);
        Gxx(n,m)=(3*(K(m,1)-P(n,1))^2/((K(m,1)-P(n,1))^2+(K(m,2)-P(n,2))^2+(K(m,3)-P(n,3))^2)-1)*(a*Dx*Dy*Dz)/((K(m,1)-P(n,1))^2+(K(m,2)-P(n,2))^2+(K(m,3)-P(n,3))^2)^1.5;
        Gxy(n,m)=(3*(K(m,1)-P(n,1))*(K(m,2)-P(n,2))/((K(m,1)-P(n,1))^2+(K(m,2)-P(n,2))^2+(K(m,3)-P(n,3))^2))*(a*Dx*Dy*Dz)/((K(m,1)-P(n,1))^2+(K(m,2)-P(n,2))^2+(K(m,3)-P(n,3))^2)^1.5;
        Gxz(n,m)=(3*(K(m,1)-P(n,1))*(K(m,3)-P(n,3))/((K(m,1)-P(n,1))^2+(K(m,2)-P(n,2))^2+(K(m,3)-P(n,3))^2))*(a*Dx*Dy*Dz)/((K(m,1)-P(n,1))^2+(K(m,2)-P(n,2))^2+(K(m,3)-P(n,3))^2)^1.5;
        Gyy(n,m)=(3*(K(m,2)-P(n,2))^2/((K(m,1)-P(n,1))^2+(K(m,2)-P(n,2))^2+(K(m,3)-P(n,3))^2)-1)*(a*Dx*Dy*Dz)/((K(m,1)-P(n,1))^2+(K(m,2)-P(n,2))^2+(K(m,3)-P(n,3))^2)^1.5;
        Gyz(n,m)=(3*(K(m,2)-P(n,2))*(K(m,3)-P(n,3))/((K(m,1)-P(n,1))^2+(K(m,2)-P(n,2))^2+(K(m,3)-P(n,3))^2))*(a*Dx*Dy*Dz)/((K(m,1)-P(n,1))^2+(K(m,2)-P(n,2))^2+(K(m,3)-P(n,3))^2)^1.5;
        Gzz(n,m)=(3*(K(m,3)-P(n,3))^2/((K(m,1)-P(n,1))^2+(K(m,2)-P(n,2))^2+(K(m,3)-P(n,3))^2)-1)*(a*Dx*Dy*Dz)/((K(m,1)-P(n,1))^2+(K(m,2)-P(n,2))^2+(K(m,3)-P(n,3))^2)^1.5;
    end
end
    
      
      