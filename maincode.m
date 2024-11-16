close all
clear all
clc

%%%%%%%%%%%%%%%%%%%%%
%Variable Definition%
%%%%%%%%%%%%%%%%%%%%%
p=10;
L=10^2;
nrx=p+1;
N_vec=[400,576,784,1024,1600];

num_realisation_H=100;
K_vec=[2,4,8,16,32,64,128];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MSE vs N  Closed Form          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
K=16;
mse_n=zeros(1,length(N_vec));
for ii=1:length(N_vec)
    N=N_vec(ii);
    ntx=1/N;
    mse_mean=0;
    M=N;
    for k=1:num_realisation_H
        H=channel(N,M,L);
        [U,S,W]=svd(H);
        W=W(:,1:K);
        [V_star,mse]=combiner(W,K,ntx,nrx,H);
        mse_mean=mse_mean+mse;
    end
mse_n(ii)=real(mse_mean/num_realisation_H);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MSE vs K                      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N_k=[400,1024];           % only consider these values of N
mse_k=zeros(length(N_k),length(K_vec));
for n=1:length(N_k)
    for m=1:length(K_vec)
        N=N_k(n);
        K=K_vec(m);
        ntx=1/N;
        mse_mean=0;
        M=N;
        for k=1:num_realisation_H
            H=channel(N,M,L);
            [U,S,W]=svd(H);
            W=W(:,1:K);
            [V_star,mse]=combiner(W,K,ntx,nrx,H);
            mse_mean=mse_mean+mse;
        end
    mse_k(n,m)=real(mse_mean/num_realisation_H);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%MSE vs N  Monte Carlo          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
K=16;
mc_realisations=2*10^3;
mse_mc_doubly=zeros(1,length(N_vec));
mse_mc_1bit=zeros(1,length(N_vec));
for ii=1:length(N_vec)
    N=N_vec(ii)
    ntx=1/N;
    mse_mean_doubly=0;
    mse_mean_1bit=0;
    M=N;
    for k=1:num_realisation_H
        H=channel(N,M,L);
        [U,S,W]=svd(H);
        W=W(:,1:K);

        %signal matrix,precoding, and Tx quantization
        S=(1/sqrt(2))*(randn(K,mc_realisations)+1i*randn(K,mc_realisations));
        X=W*S;
        T=sqrt(ntx/2)*(sign(real(X))+1i*sign(imag(X)));
        
        %Nose matrix
        Z=(1/sqrt(2))*(randn(M,mc_realisations)+1i*randn(M,mc_realisations));
        
        %Rx matrix and processing for 1-bit ADC and doubly ADC/DAC
        Y_doubly=sqrt(p)*H*T+Z;
        Y_1bit=sqrt(p)*H*X+Z;
        R_doubly=sqrt(nrx/2)*(sign(real(Y_doubly))+1i*sign(imag(Y_doubly)));
        R_1bit=sqrt(nrx/2)*(sign(real(Y_1bit))+1i*sign(imag(Y_1bit)));

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Estimate Cr and E[rs'] for Combiner Design%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        Crr_doubly=zeros(M,M);
        Crs_doubly=zeros(M,K);
        Crr_1bit=zeros(M,M);
        Crs_1bit=zeros(M,K);
        for n=1:mc_realisations
            Crr_doubly=Crr_doubly+R_doubly(:,n)*R_doubly(:,n)';
            Crr_1bit=Crr_1bit+R_1bit(:,n)*R_1bit(:,n)';
            Crs_doubly=Crs_doubly+R_doubly(:,n)*S(:,n)';
            Crs_1bit=Crs_1bit+R_1bit(:,n)*S(:,n)';
        end
        Crr_doubly=Crr_doubly/mc_realisations;
        Crs_doubly=Crs_doubly/mc_realisations;
        Crs_1bit=Crs_1bit/mc_realisations;
        Crr_1bit=Crr_1bit/mc_realisations;

        %Combiner computation
        V_star_doubly=inv(Crr_doubly)*Crs_doubly;
        V_star_1bit=inv(Crr_1bit)*Crs_1bit;
        S_hat_doubly=V_star_doubly'*R_doubly;
        S_hat_1bit=V_star_1bit'*R_1bit;

        %Error Processing
        E_vec_doubly=S-S_hat_doubly;
        E_vec_doubly=sum(abs(E_vec_doubly).^2);
        E_vec_1bit=S-S_hat_1bit;
        E_vec_1bit=sum(abs(E_vec_1bit).^2);
        mse_mean_doubly=mse_mean_doubly+(1/K)*mean(E_vec_doubly);
        mse_mean_1bit=mse_mean_1bit+(1/K)*mean(E_vec_1bit);
    end
    mse_mc_doubly(ii)=real(mse_mean_doubly/num_realisation_H);
    mse_mc_1bit(ii)=real(mse_mean_1bit/num_realisation_H);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Data Detections for 16-PSK data%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N_Data=[400,1024,1600];        %only consider these values of N
K=8;
data_vec=[];
decode_vec_doubly=[];
sym_sent_vec=[];
sym_recv_vec_doubly=[];
for n=1:length(N_Data)
    N=N_Data(n);
    M=N;
    H=channel(N,M,L);
    [U,S,W]=svd(H);
    W=W(:,1:K);
   [V_star,mse]=combiner(W,K,ntx,nrx,H);

    data_vec=[];
    decode_vec_doubly=[];
    sym_sent_vec=[];
    sym_recv_vec_doubly=[];
    sym_recv_vec_1bit=[];
    decode_vec_1bit=[];
    figure(n*100)
    %Loop through 100000 data realization for ser estimation
    for k=1:100000

        %16 PSK Modulation
        data=(randi([0,15],K,1));
        s = pskmod(data,16);
        sym_sent_vec=[sym_sent_vec;s];
        data_vec=[data_vec; data];
        
        %Tx Processing 
        x=W*s;
        t=sqrt(ntx/2)*(sign(real(x))+1i*sign(imag(x)));

        %AWGN 
        z=(1/sqrt(2))*(randn(M,1)+1i*randn(M,1));

        %Rx processing for 1-bit ADC and doubly ADC/DAC
        y_doubly=sqrt(p)*H*t+z;
        y_1bit=sqrt(p)*H*x+z;
        r_doubly=sqrt(nrx/2)*(sign(real(y_doubly))+1i*sign(imag(y_doubly)));
        r_1bit=sqrt(nrx/2)*(sign(real(y_1bit))+1i*sign(imag(y_1bit)));
        soft_doubly=V_star'*r_doubly;
        soft_1bit=V_star'*r_1bit;
        sym_recv_vec_1bit=[sym_recv_vec_1bit;soft_1bit];
        sym_recv_vec_doubly=[sym_recv_vec_doubly;soft_doubly];

        %Demodulation and Decode
        hard_doubly=pskdemod(soft_doubly,16);
        hard_1bit=pskdemod(soft_1bit,16);
        decode_vec_doubly=[decode_vec_doubly; hard_doubly];
        decode_vec_1bit=[decode_vec_1bit; hard_1bit];
    end
    
    %%%%calculate ser
    parity_doubly=decode_vec_doubly~=data_vec;
    parity_1bit=decode_vec_1bit~=data_vec;
    ser_doubly=sum(parity_doubly)/length(data_vec);
    ser_1bit=sum(parity_1bit)/length(data_vec);

    %constellation
    ind=randi([1,100000],2000,1);
    hold on
    plot(real(sym_recv_vec_doubly(ind)),imag(sym_recv_vec_doubly(ind)),'ko')
    hold on
    plot(real(sym_recv_vec_1bit(ind)),imag(sym_recv_vec_1bit(ind)),'bo')
    hold on
    plot(real(sym_sent_vec(ind)),imag(sym_sent_vec(ind)),'ro',...
        'LineWidth',2)

    ylabel("Q")
    xlabel("I")
    hold off


set(legend( strcat('$$\hat{s}/$$(1-bit ADCs/DACs): SER = ',...
sprintf('%10e',ser_doubly)),strcat('$$\hat{s}/$$(1-bit ADCs): SER = ',...
sprintf('%10e',ser_1bit)),'s(sent)'),'Interpreter','Latex','FontSize', 10)

xlabel('I')
ylabel('Q')
title(strcat("N = M = ",num2str(N), ", K=8, p =10 "))
end


%%%%%%%%%%%%%%%%%
plot generations%
%%%%%%%%%%%%%%%%%
figure(1)
plot(N_vec,mse_n)
hold on
plot(N_vec,mse_mc_doubly,"*-")
hold on
plot(N_vec,mse_mc_1bit,"-.")
legend('MSE approx (1bit-DAC/ADC)','MSE Monte Carlo (1bit-DAC/ADC)',...
    'MSE Monte Carlo (1bit-ADC')
title('K = 16, p = 10 dB')
xlabel("N = M")
ylabel("MSE")


figure(2)
plot(K_vec,mse_k(1,:),"--o")
hold on
plot(K_vec,mse_k(2,:),"-.*")
legend('N=400','N=1024')
title("p = 10 dB")
xlabel("K")
ylabel("MSE approx")
xticks(0:8:128)
hold off

