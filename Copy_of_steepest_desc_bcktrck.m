function [delta, fk] = Copy_of_steepest_desc_bcktrck(x0, f, gradf, alpha0, kmax, tolgrad, c1, rho, btmax,xmin)
    xk = x0;
    fk = f(xk);
    for kk = 1:kmax
        pk = -gradf(xk); %descent direction
        check = true;
        alpha = alpha0;
        t = 0;
        xold=xk;
        % backtracking
        while check
            if t == btmax-1
                check = false;
            end
            xk = xold+alpha*pk;
            temp = f(xk);
            if temp <= fk-c1*alpha*pk'*pk
                fk = temp;
                check = false;
            else
                alpha = rho*alpha;
                t = t+1;
            end
        end
        gradfk_norm = norm(pk);
        if gradfk_norm < tolgrad
            break;
        end
    end
    delta=norm(xk-xmin,2);
end