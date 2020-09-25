function[E] = Calc_EA(M,e)

    tol = 10^-8;
    Etemp = M;
    ratio = 1;

    while abs(ratio) > tol
        f_E = Etemp - e*sin(Etemp) - M;
        f_Eprime = 1 - e*cos(Etemp);
        ratio = f_E/f_Eprime;
        if abs(ratio) > tol
            Etemp = Etemp - ratio;
        else
            E = Etemp;
        end
    end
    
end
