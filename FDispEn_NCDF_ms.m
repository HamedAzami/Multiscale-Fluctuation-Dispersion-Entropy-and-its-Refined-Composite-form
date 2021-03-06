function [Out_FDispEn npdf]=FDispEn_NCDF_ms(x,m,nc,mu,sigma,tau)
%
% This function calculates fluctuation-based dispersion entropy (FDispEn) of a univariate
% signal x, using normal cumulative distribution function (NCDF) with defined mean (mu) and 
% standard deviation (sigma) values.
%
% Inputs:
%
% x: univariate signal - a vector of size 1 x N (the number of sample points)
% m: embedding dimension
% nc: number of classes (it is usually equal to a number between 3 and 9 - we used c=6 in Ref. [1])
%mu: mean value
%sigma: standard deviation. 
% tau: time lag (it is usually equal to 1)
%
% Outputs:
%
% Out_DispEn: scalar quantity - the DispEn of x
% npdf: a vector of length (2*nc-1)^(m-1), showing the normalized number of fluctuation-based disersion patterns of x
%
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

N=length(x);


%% Mapping approaches

y=normcdf(x,mu,sigma);
y(y==1)=1-1e-10;
y(y==0)=1e-10;
z=round(y*nc+0.5);

all_patterns=[1:2*nc-1]';

for f=2:m-1
    temp=all_patterns;
    all_patterns=[];
    j=1;
    for w=1:2*nc-1
        [a,b]=size(temp);
        all_patterns(j:j+a-1,:)=[temp,w*ones(a,1)];
        j=j+a;
    end
end


N_PDE=(2*nc-1)^(m-1);
for i=1:N_PDE
    key(i)=0;
    for ii=1:m-1
        key(i)=key(i)*100+all_patterns(i,ii);
    end
end


for i_h=1:m
    ind(i_h,:)=[(i_h-1)*tau+1:N-(m-1)*tau+(i_h-1)*tau];
end

embd2 = z(ind(:,1:end));
dembd2=diff(embd2)'+nc;

emb=zeros(N-(m-1)*tau,1);
for i_e=m-1:-1:1
    emb=dembd2(:,i_e)*100^(i_e-1)+emb;
end

pdf=zeros(1,N_PDE);

for id=1:N_PDE
    [R,C]=find(emb==key(id));
    pdf(id)=length(R);
end

npdf=pdf/(N-(m-1)*tau);
p=npdf(npdf~=0);
Out_FDispEn = -sum(p .* log(p));