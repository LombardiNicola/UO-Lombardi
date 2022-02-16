function [xk, fk, gradfk_norm, k, xseq, btseq] = ...
    newton_bcktrck(x0, f, gradf, Hessf, kmax, tolgrad, c1, rho, btmax)
    xk = x0;
    fk = f(xk);
    xseq = zeros(length(x0), kmax+1); % steps
    btseq = zeros(1, kmax+1); % number of backtracking iterations
    xseq(:, 1) = xk;
    for k = 1:kmax
        gradfk = gradf(xk);
        Hessfk = Hessf(xk);
        pk = Hessfk\-gradfk; % descent direction
        alpha = 1;
        check = true; % stopping inner cycle
        t = 0; % inner cycle iterator
        % backtraking
        while check
            if t == btmax-1
                check = false;
            end
            xk = xseq(:, k) + alpha*pk;
            temp = f(xk);
            if temp <= fk - c1*alpha*pk'*pk
                check = false;
            else
                alpha = rho*alpha;
                t = t+1;
            end
        end
        fk = temp;
        xseq(:, k+1) = xk;
        btseq(k+1) = t;
        gradfk_norm = norm(pk);
        % exit condition
        if gradfk_norm < tolgrad
            break;
        end
    end
    xseq = xseq(:, 1:k+1);
    btseq = btseq(1:k+1);
end