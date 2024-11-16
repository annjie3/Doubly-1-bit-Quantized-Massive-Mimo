function [V_star,mse] = combiner(W,K,ntx,nrx,H)
%This funtion generates the combiner needed for data detections
p=10;
K=size(W,2);
N=size(W,1);
M=N;

Cx=W*W';

%%Compute Ct
diagCx=diag(Cx);
diagCx_pow_nhalf=1./sqrt(diagCx);
Cx_nhalf_mat=diag(diagCx_pow_nhalf);

Ct=(2/pi)*ntx*(asin(Cx_nhalf_mat*real(Cx)*Cx_nhalf_mat)+...
    1i*asin(Cx_nhalf_mat*imag(Cx)*Cx_nhalf_mat));


%%%Compute Cy
Cy=p*H*Ct*H'+eye(M);
diagCy=diag(Cy);
diagCy_pow_nhalf=1./sqrt(diagCy);
Cy_nhalf_mat=diag(diagCy_pow_nhalf);

Cr=(2/pi)*nrx*(asin(Cy_nhalf_mat*real(Cy)*Cy_nhalf_mat)+...
    1i*asin(Cy_nhalf_mat*imag(Cy)*Cy_nhalf_mat));

%%%Compute Gtx and Grx
Gtx=sqrt(2*ntx/pi)*Cx_nhalf_mat;
Grx=sqrt(2*nrx/pi)*Cy_nhalf_mat;

%Compute V_star
V_star=sqrt(p)*inv(Cr+(1e-14)*eye(length(Cr)))*Grx*H*Gtx*W;
mse=1+(1/K)*trace(V_star'*Cr*V_star)-...
    (2/K)*sqrt(p)*trace(real(V_star'*Grx*H*Gtx*W));
end