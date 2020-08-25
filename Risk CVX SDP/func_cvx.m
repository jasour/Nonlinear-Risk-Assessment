function [XX]=func_cvx(A,C,b,n,m)

cvx_begin
    variable XX( n, n ) symmetric
    dual variables y{m}
sdpsettings('solver','mosek');

    minimize( trace(C*XX) );
       for k = 1 : m;
        trace(A{k}*XX)==b(k) : y{k};  % coeff of x^k , b=[x^0,x,...]
    end
    XX == semidefinite(n);
cvx_end

end
