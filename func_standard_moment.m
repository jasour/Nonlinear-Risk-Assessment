function [Ap,Mind]=func_standard_moment(nvar,d)

% To construct Ap, we use the moment matirx structure.
% Lecture 5: Duality of SOS and Moment based Semidefinite Programs (SDPs),page 125
% https://rarnop.mit.edu/Lectures-Codes

Ny= round(factorial(nvar+2*d)/(factorial(2*d)*factorial(nvar)));   
vpow=[]; for k = 0:d; vpow = [vpow;genpow(nvar,k)]; end


B=sparse(size(vpow,1),size(vpow,1));
Ap=cell(1,Ny); Ap(:) = {B};

k=0;
for i=1:size(vpow,1)   
    for j=1:i        
       k=k+1;clc;disp('Moment MAptrix'); disp([i,j,size(vpow,1)])    
       Mind(i,j)=(glex2num(vpow(i,:)+vpow(j,:)));Mind(j,i)=Mind(i,j);
       Ap{glex2num(vpow(i,:)+vpow(j,:))}(i,j)=1;
       Ap{glex2num(vpow(i,:)+vpow(j,:))}(j,i)=1;
    end
end



end