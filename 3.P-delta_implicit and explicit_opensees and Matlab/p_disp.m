clc;clear all;close all;
A0=5;L0=160;E=36000000;N0=1500;R=200;Kel=A0*E/L0;


%a_1=[-cos(teta_1),-sin(teta_1),cos(teta_1),sin(teta_1)]
%a_2=[-cos(teta_2),-sin(teta_2),cos(teta_2),sin(teta_2)]

%teta_2=-teta_1
%sin(teta_1)=u/L
%v1=u*a_1(1,4)
%pint_1=N*a_1(1,4)
%sin(teta_2)=-sin(teta_1)
%v2=-u*a_2(1,2)
%pint_2=N*a_1(1,2)


%k1=[(cos(teta_1))^2,cos(teta_1)*sin(teta_1);(sin(teta_1))^2,cos(teta_1)*sin(teta_1)]
%k2=[(cos(teta_2))^2,cos(teta_2)*sin(teta_2);(sin(teta_2))^2,cos(teta_2)*sin(teta_2)]

%KL_1=Kel*[k1,-k1;-k1,k1]
%KL_2=Kel*[k2,-k2;-k2,k2]
 %K_T1=KL_1(4,4)+Knl(4,4);
     %   K_T2=KL_2(2,2)+Kn1(2,2);
     %   K_T(i)=K_T1+K_T2; 
%j=[1,1,1,1];
%f=[-1,-1];
%Knl=(N/L)*((diag(j)+diag(f,2)+diag(f,-2));

%pint=pint_1+pint_2

i=0;p_int=0;N=N0;L=L0;sin=0;u=0;v=0;


tolerance=10^-4
step=10
for  s=1:step
    p_ext=(s/step)*R
    
    p_r=p_ext-p_int(i+1)
    
    while  p_r>tolerance
          
        i=i+1
        
        KT(i)=(2*A0*E/L0)*(sin^2)+(2*N/L)
        
        delta_u=p_r/KT(i);
        
        u(i+1)=u(i)+delta_u;
        
        L=sqrt((u(i+1)^2)+(L0^2));
        
        sin=u(i+1)/L;
        
        delta_v=sin*delta_u
        
        v(i+1)=v(i)+delta_v
        
        p_int(i+1)=2*sin*(N0+(A0*E/L0)*v(i+1))
        
       p_r=p_ext-p_int(i+1)
       
        % A(1:2,i+1)=p_int'
       %  B(1:2,i+1)=u'
        plot(u,p_int);grid on; hold on
    end
end


        
        
    
    
    




