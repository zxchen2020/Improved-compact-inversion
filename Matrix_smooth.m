function W_xyz=Matrix_smooth(M,Mx,My,Mz)
WZ=sparse(M,M);WX=sparse(M,M);WY=sparse(M,M);
for k=1:My 
for j=1:Mx 
for i=1:(Mz-1)  
    WZ((k-1)*Mx*Mz+(j-1)*Mz+i,(k-1)*Mx*Mz+(j-1)*Mz+i)=-1;
    WZ((k-1)*Mx*Mz+(j-1)*Mz+i,(k-1)*Mx*Mz+(j-1)*Mz+i+1)=1;
end
end
end
for k=1:My 
for j=1:(Mx-1)
for i=1:Mz  
    WX((k-1)*Mx*Mz+(j-1)*Mz+i,(k-1)*Mx*Mz+(j-1)*Mz+i)=-1;
    WX((k-1)*Mx*Mz+(j-1)*Mz+i,(k-1)*Mx*Mz+(j-1)*Mz+i+Mz)=1;  
end
end
end
for k=1:(My-1) 
for j=1:Mx
for i=1:Mz  
    WY((k-1)*Mx*Mz+(j-1)*Mz+i,(k-1)*Mx*Mz+(j-1)*Mz+i)=-1;
    WY((k-1)*Mx*Mz+(j-1)*Mz+i,(k-1)*Mx*Mz+(j-1)*Mz+i+Mz*Mx)=1;  
end
end
end
W_xyz=WX'*WX+WY'*WY+WZ'*WZ;%¹â»¬Ô¼Êø¾ØÕó

