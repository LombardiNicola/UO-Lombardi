function [delta, fk] = ...
    Copy_of_newton_bcktrck(x0, f, gradf, Hessf, kmax, tolgrad, c1, rho, btmax,xmin)
    xk = x0;
    fk = f(xk);
    for k = 1:kmax
        gradfk = gradf(xk);
        Hessfk = Hessf(xk);
        pk = Hessfk\-gradfk; % descent direction
        alpha = 1;
        xold=xk;
        xk = xk + alpha*pk;
        check = true;
        t = 0;
        % backtraking
        while check
            if t == btmax-1
                check = false;
            end
            xk = xold + alpha*pk;
            temp = f(xk);
            if temp <= fk - c1*alpha*pk'*pk
                check = false;
            else
                alpha = rho*alpha;
                t = t+1;
            end
        end
        fk = temp;
        gradfk_norm = norm(pk);
        if gradfk_norm < tolgrad
            break;
        end
    end
    delta=norm(xk-xmin,2);
end