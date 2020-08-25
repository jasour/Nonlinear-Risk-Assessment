function [As1g]=func_standard_localizing(g,nvar,d)

% To construct As1g, we use the localizing matirx structure.
% Lecture 5: Duality of SOS and Moment based Semidefinite Programs (SDPs), page 129
% https://rarnop.mit.edu/Lectures-Codes
clc
Ny= round(factorial(nvar+2*d)/(factorial(2*d)*factorial(nvar))); % Number of y variables : Moments of Measure Mu  

Cg=g.coef; Cc=0+Cg(:,:); % Cc1 : coefficents of g(x)
Dg=g.deg;  Dc=0+Dg(:,:); % Bc1 : degrees of monomials of g(x)
d1=floor(abs(2*d- max(sum(Dc,2)) )/2);


n2=size(Dc,2);
vpow=[]; for k = 0:d1; vpow = [vpow;genpow(n2,k)]; end; size(vpow,1)

B=sparse(size(vpow,1),size(vpow,1));
As1g=cell(1,Ny); As1g(:) = {B};


for i=1:size(vpow,1)   
    for j=1:i
        clc;disp('Localization Matrix');disp([i,j,size(vpow,1)])
        a=vpow(i,:)+vpow(j,:);   
        for k=1:size(Dc)                        
            B=As1g{glex2num(a+Dc(k,:))};
            B(i,j)=Cc(k);B(j,i)=Cc(k);
            As1g{glex2num(a+Dc(k,:))}=B;
            As1g{glex2num(a+Dc(k,:))}(i,j)=Cc(k);
            As1g{glex2num(a+Dc(k,:))}(j,i)=Cc(k);                 
        end        
    end
end

end
