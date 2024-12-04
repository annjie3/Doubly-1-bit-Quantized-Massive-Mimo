function [H] = channel(N,M,L)
%This function genrates a channel based on discrete physical model use 
%in the paper based on [15] for UPA arrays

scatter_dim=pi/12;
theta_vec_T=scatter_dim*(-1+2*randn(L,1));
phi_vec_T=scatter_dim*(-1+2*randn(L,1));

theta_vec_R=scatter_dim*(-1+2*randn(L,1));
phi_vec_R=scatter_dim*(-1+2*randn(L,1));
beta=(1/sqrt(2))*(randn(L,1)+1i*randn(L,1));

Hp= diag(N*beta/sqrt(L));
AT=zeros(N,L);
AR=zeros(M,L);

for m=1:L
    theta_T=theta_vec_T(m);
    phi_T=phi_vec_T(m);
    theta_R=theta_vec_R(m);
    phi_R=phi_vec_R(m);
    
    bt=sin(theta_T)*cos(phi_T)*[0:sqrt(N)-1]'+...
        sin(theta_T)*sin(phi_T)*[0:sqrt(N)-1];

    br=sin(theta_R)*cos(phi_R)*[0:sqrt(M)-1]'+...
        sin(theta_R)*sin(phi_R)*[0:sqrt(M)-1];

    aT=exp(1i*pi*reshape(bt',[],1));
    aR=exp(1i*pi*reshape(br',[],1));
    AT(:,m)=aT;
    AR(:,m)=aR;
end

H=AR*Hp*AT';
end