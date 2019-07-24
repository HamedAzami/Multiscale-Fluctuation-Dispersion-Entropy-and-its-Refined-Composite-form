function Out_RCMFDE=RCMFDE(x,m,c,tau,Scale)
%
% This function calculates the refined composite multiscale fluctuation-based dispersion entropy (RCMFDE) of a univariate signal x
%
% Inputs:
%
% x: univariate signal - a vector of size 1 x N (the number of sample points)
% m: embedding dimension
% nc: number of classes (it is usually equal to a number between 3 and 9 - we used c=6 in Ref. [1])
% tau: time lag (it is usually equal to 1)
% Scale: maximum scale factor
%
% Outputs:
%
% Out_RCMFDE: scalar quantity - the RCMFDE of x
%
%
% EXAMPLE: RCMFDE(rand(1,1000),2,6,1,10)
%
% Ref:
% [1] H. Azami, S. Arnold, S. Sanei, Z. Chang, G. Sapiro, J. Escudero, and A. Gupta, "Multiscale Fluctuation-based
% Dispersion Entropy and its Applications to Neurological Diseases", IEEE ACCESS, 2019.
% [2] H. Azami, and J. Escudero, "Amplitude-and Fluctuation-Based Dispersion Entropy", Entropy, vol. 20, no. 3, p.210, 2018.
%
% If you use the code, please make sure that you cite the references [1] and [2].
%
% Hamed Azami
% hazami@mgh.harvard.edu and hmd.azami@gmail.com
%
%  17-April-2019
%%

Out_RCMFDE=NaN*ones(1,Scale);

Out_RCMFDE(1)=FDispEn_NCDF(x,m,c,tau);

sigma=std(x);
mu=mean(x);

for j=2:Scale
    pdf=[];
    for jj=1:j
        xs = Multi(x(jj:end),j);
        [DE, T_pdf]=FDispEn_NCDF_ms(xs,m,c,mu,sigma,tau);
        pdf=[pdf ; T_pdf];
    end
    pdf=mean(pdf,1);
    pdf=pdf(pdf~=0);
    Out_RCMFDE(j)=-sum(pdf .* log(pdf));
end


function M_Data = Multi(Data,S)

%  generate the consecutive coarse-grained time series
%  Input:   Data: time series;
%           S: the scale factor
% Output:
%           M_Data: the coarse-grained time series at the scale factor S

L = length(Data);
J = fix(L/S);

for i=1:J
    M_Data(i) = mean(Data((i-1)*S+1:i*S));
end
