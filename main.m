%% Rosenbrock steepest descent%%
fRosenbrock = @(x) 100*(x(2)-x(1)^2)^2 + (1-x(1))^2;
gradRos = @(x) [-400*(x(2)-x(1)^2)*x(1)-2*(1-x(1)); 200*(x(2)-x(1)^2)];
rho=0.5;
alpha0=1;
kmax=10000;
tolgrad=1e-8;
btmax=50;
c=1e-4;
f=fRosenbrock;
gradf=gradRos;

x0_A = [1.2 1.2];
x0_B = [-1.2 1];

x0=x0_A;
[xk_A, fk_A, gradfk_norm_A, kk_A, xseq_A, btseq_A] = ...
    steepest_desc_bcktrck(x0, f, gradf, alpha0, kmax, tolgrad, c, rho, btmax);
x0=x0_B;
[xk_B, fk_B, gradfk_norm_B, kk_B, xseq_B, btseq_B] = ...
    steepest_desc_bcktrck(x0, f, gradf, alpha0, kmax, tolgrad, c, rho, btmax);

figure('Position', [100 100 1000 330])
red = [1 0 0]; % darker color are used for functions
orange = [1 0.6  0]; % lighter for the gradient norm
blue = [0 0 1]; 
l_blue = [0.7 0.7 1];
hold on
yyaxis left
set(gca,'YColor', [0 0 0]);
fseq_A = arrayfun(@(i) f(xseq_A(:,i)), 1:kk_A+1);
plot(1:kk_A+1, fseq_A, '-', "Color", red, "LineWidth", 2);
fseq_B = arrayfun(@(i) f(xseq_B(:,i)), 1:kk_B+1);
plot(1:kk_B+1, fseq_B, '-', "Color", blue, "LineWidth", 2);
yyaxis right
set(gca,'YColor', [0 0 0]);
gradNormseq_A = arrayfun(@(i) norm(gradf(xseq_A(:,i)),2), 1:kk_A+1);
plot(1:kk_A+1, gradNormseq_A, '-', "Color", orange, "LineWidth", 2);
gradNormseq_B = arrayfun(@(i) norm(gradf(xseq_B(:,i)),2), 1:kk_B+1);
plot(1:kk_B+1, gradNormseq_B, '-', "Color", l_blue, "LineWidth", 2);
legend ({'f case A', 'f case B', ...
    'gradf case A', 'gradf case B'}) 
hold off
saveas(gcf, 'fGradfSteepRos.png')

figure('Position', [10 10 1000 1000])
fsurf(@(x,y) f([x y]), [-2 2 -2 2],'FaceAlpha',0.4);
hold on
set(gca,'YColor', [0 0 0]);
plot3(xseq_A(1,:),xseq_A(2,:),fseq_A, '-','LineWidth', 2,'Color', red);
plot3(xseq_B(1,:),xseq_B(2,:),fseq_B, '-','LineWidth', 2,'Color', blue);
hold off
% saveas(gcf, 'surfSteepRos.png')
%% Rosenbrock Newton%%
fRosenbrock = @(x) 100* ( x(2)-x(1)^2 )^2 + (1-x(1))^2;
gradRos = @(x) [-400*(x(2)-x(1)^2)*x(1)-2*(1-x(1)); 200*(x(2)-x(1)^2)];
HessRos = @(x) [ -400*(-2*x(1)^2+x(2)-x(1)^2)+2 -400*x(1) ; -400*x(1) 200];
rho=0.5;
kmax=10000;
tolgrad=1e-8;
btmax=100;
c=1e-4;
f=fRosenbrock;
gradf=gradRos;
Hessf=HessRos;
x0_A = [1.2 1.2];
x0_B = [-1.2 1];

x0=x0_A;
[xk_A, fk_A, gradfk_norm_A, kk_A, xseq_A, btseq_A] = ...
    newton_bcktrck(x0, f, gradf, Hessf, kmax, tolgrad, c, rho, btmax);
x0=x0_B;
[xk_B, fk_B, gradfk_norm_B, kk_B, xseq_B, btseq_B] = ...
    newton_bcktrck(x0, f, gradf, Hessf, kmax, tolgrad, c, rho, btmax);

figure('Position', [100 100 1000 330])
red = [1 0 0]; % darker color are used for functions
orange = [1 0.6  0]; % lighter for the gradient norm
blue = [0 0 1]; 
l_blue = [0.7 0.7 1];
hold on
yyaxis left
set(gca,'YColor', [0 0 0]);
fseq_A = arrayfun(@(i) f(xseq_A(:,i)), 1:kk_A+1);
plot(1:kk_A+1, fseq_A, '-', "Color", red, "LineWidth", 2);
fseq_B = arrayfun(@(i) f(xseq_B(:,i)), 1:kk_B+1);
plot(1:kk_B+1, fseq_B, '-', "Color", blue, "LineWidth", 2);
set(gca,'YColor', [0 0 0]);
yyaxis right
set(gca,'YColor', [0 0 0]);
gradNormseq_A = arrayfun(@(i) norm(gradf(xseq_A(:,i)),2), 1:kk_A+1);
plot(1:kk_A+1, gradNormseq_A, '-', "Color", orange, "LineWidth", 2);
gradNormseq_B = arrayfun(@(i) norm(gradf(xseq_B(:,i)),2), 1:kk_B+1);
plot(1:kk_B+1, gradNormseq_B, '-', "Color", l_blue, "LineWidth", 2);
legend ({'f case A', 'f case B', ...
    'gradf case A', 'gradf case B'}) 
set(gca,'YColor', [0 0 0]);
hold off
saveas(gcf, 'fGradfNewtonRos.png')

figure('Position', [10 10 1000 1000])
fsurf(@(x,y) f([x y]), [-2 2 -2 2],'FaceAlpha',0.4);
hold on
set(gca,'YColor', [0 0 0]);
plot3(xseq_A(1,:),xseq_A(2,:),fseq_A, '-','LineWidth', 2,'Color', red);
plot3(xseq_B(1,:),xseq_B(2,:),fseq_B, '-','LineWidth', 2,'Color', blue);
hold off
% saveas(gcf, '.png')
%% Chained Rosenbrock steepest %%
rho=0.5;
alpha0=1;
kmax=2000;
tolgrad=1e-8;
btmax=50;
c=1e-4;
f=@(x) fCR(x);
gradf=@(x) gCR(x);


ns=[10:10:90 100:100:1000];
len=length(ns);
rounds=10;
times=zeros(len,2);
delta_ij=zeros(len,rounds);
f_ij=zeros(len,rounds);
delta2_ij=zeros(len,rounds);
f2_ij=zeros(len,rounds);
for i=1:len
    i
    n=ns(i);
    x0 = 1.2*ones(n,1);
    xmin=ones(n,1);
    for j=1:rounds
        t=tolgrad*sqrt(n);
        tic
        [delta_x, fk] = Copy_of_steepest_desc_bcktrck(x0, f, gradf, alpha0, kmax, t, c, rho, btmax,xmin);
        times(i,1)=(j-1)/j*times(i,1)+1/j*toc;
        delta_ij(i,j)=delta_x;
        f_ij(i,j)=fk;
    end
    x0=arrayfun(@(i) 1-2.2*mod(i,2),1:n)';
    for j=1:rounds
        t=tolgrad*sqrt(n);
        tic
        [delta_x, fk] = Copy_of_steepest_desc_bcktrck(x0, f, gradf, alpha0, kmax, t, c, rho, btmax,xmin);
        times(i,2)=(j-1)/j*times(i,2)+1/j*toc;
        delta2_ij(i,j)=delta_x;
        f2_ij(i,j)=fk;
    end
end
deltas=sum(delta_ij,2)/rounds;
deltas2=sum(delta2_ij,2)/rounds;



%% Graphs %%
figure('Position', [100 100 1000 330])
hold on
yyaxis left
set(gca,'YColor', red);
plot(ns,times(:,1), '-','LineWidth', 2,'Color', red);
yyaxis right
set(gca,'YColor', blue);
plot(ns,times(:,2), '-','LineWidth', 2,'Color', blue);
legend ({'elapsed time case A', 'elapsed time case B'}) 
hold off
saveas(gcf, 'timeCRsteep.png')

figure('Position', [100 100 1000 330])
hold on

yyaxis left
plot(ns,deltas, '-','LineWidth', 2,'Color', red)
plot(ns,deltas2, '-','LineWidth', 2,'Color', blue)
set(gca,'YColor', [0 0 0]);
yyaxis right
plot(ns,arrayfun(@(i) deltas(i)/ns(i),1:19), '-','LineWidth', 2,'Color', orange)
plot(ns,arrayfun(@(i) deltas2(i)/ns(i),1:19), '-','LineWidth', 2,'Color', l_blue)
legend ({'distance from minimum case A', 'distance from minimum case B', ...
    'distance from minimum/components case A', 'distance from minimum/components case B'}) 
set(gca,'YColor', [0 0 0]);
hold off
saveas(gcf, 'distCRsteep.png')

%% Chained Rosenbrock newton %%
fRosenbrock = @(x) 100*(x(2)-x(1)^2)^2 + (1-x(1))^2;
gradRos = @(x) [-400*(x(2)-x(1)^2)*x(1)-2*(1-x(1)); 200*(x(2)-x(1)^2)];
fChainRosenbrock = @(x) fCR(x);
gradChainRos = @(x) gCR(x);
hessChainRos = @(x) hCR(x);
rho=0.5;
alpha0=1;
kmax=2000;
tolgrad=1e-8;
btmax=50;
c=1e-4;
f=fChainRosenbrock;
gradf=gradChainRos;
Hessf=hessChainRos;


ns=[10:10:90 100:100:1000];
len=length(ns);
rounds=1;
times=zeros(len,2);
delta_ij=zeros(len,rounds);
f_ij=zeros(len,rounds);
delta2_ij=zeros(len,rounds);
f2_ij=zeros(len,rounds);
for i=1:len
    i
    n=ns(i);
    x0 = 1.2*ones(n,1);
    xmin=ones(n,1);
    for j=1:rounds
        t=tolgrad*sqrt(n);
        tic
        [delta_x, fk] = Copy_of_newton_bcktrck(x0, f, gradf, Hessf, kmax, t, c, rho, btmax,xmin);
        times(i,1)=(j-1)/j*times(i,1)+1/j*toc;
        delta_ij(i,j)=delta_x;
        f_ij(i,j)=fk;
    end
    x0=arrayfun(@(i) 1-2.2*mod(i,2),1:n)';
    for j=1:rounds
        t=tolgrad*sqrt(n);
        tic
        [delta_x, fk] = Copy_of_newton_bcktrck(x0, f, gradf, Hessf, kmax, t, c, rho, btmax,xmin);
        times(i,2)=(j-1)/j*times(i,2)+1/j*toc;
        delta2_ij(i,j)=delta_x;
        f2_ij(i,j)=fk;
    end
end
deltas=sum(delta_ij,2)/rounds;
deltas2=sum(delta2_ij,2)/rounds;



%% Graphs %%
figure('Position', [100 100 1000 330])
hold on
yyaxis left
set(gca,'YColor', red);
plot(ns,times(:,1), '-','LineWidth', 2,'Color', red);
yyaxis right
set(gca,'YColor', blue);
plot(ns,times(:,2), '-','LineWidth', 2,'Color', blue);
legend ({'elapsed time case A', 'elapsed time case B'}) 
hold off
saveas(gcf, 'timeCRnewton.png')

figure('Position', [100 100 1000 330])
hold on
yyaxis left
plot(ns,deltas, '-','LineWidth', 2,'Color', red)
plot(ns,deltas2, '-','LineWidth', 2,'Color', blue)
set(gca,'YColor', [0 0 0]);
yyaxis right
plot(ns,arrayfun(@(i) deltas(i)/ns(i),1:19), '-','LineWidth', 2,'Color', orange)
plot(ns,arrayfun(@(i) deltas2(i)/ns(i),1:19), '-','LineWidth', 2,'Color', l_blue)
legend ({'distance from minimum case A', 'distance from minimum case B', ...
    'distance from minimum/components case A', 'distance from minimum/components case B'}) 
set(gca,'YColor', [0 0 0]);
hold off
saveas(gcf, 'distCRnewton.png')

%% Brown steepest %%
f = @(x) fB(x);
gradf = @(x) gB(x);
rho=0.5;
alpha0=1;
kmax=5000;
tolgrad=1e-8;
btmax=50;
c=1e-4;


ns=[10:10:90 100:100:1000];
len=length(ns);
rounds=1;
times=zeros(len,2);
delta_ij=zeros(len,rounds);
f_ij=zeros(len,rounds);
delta2_ij=zeros(len,rounds);
f2_ij=zeros(len,rounds);
for i=1:len
    i
    n=ns(i);
    x0 = arrayfun(@(i) 1+mod(i,2)  ,1:n);
    xmin = arrayfun(@(i) 3+log(20)/20+mod(i,2)*(3-(3+log(20)/20)),1:n);
    for j=1:rounds
        t=tolgrad*sqrt(n);
        tic
        [delta_x, fk] = Copy_of_steepest_desc_bcktrck(x0, f, gradf, alpha0, kmax, t, c, rho, btmax,xmin);
        times(i,1)=(j-1)/j*times(i,1)+1/j*toc;
        delta_ij(i,j)=delta_x;
        f_ij(i,j)=fk;
    end
    x0=zeros(n,1);
    for j=1:rounds
        t=tolgrad*sqrt(n);
        tic
        [delta_x, fk] = Copy_of_steepest_desc_bcktrck(x0, f, gradf, alpha0, kmax, t, c, rho, btmax,xmin);
        times(i,2)=(j-1)/j*times(i,2)+1/j*toc;
        delta2_ij(i,j)=delta_x;
        f2_ij(i,j)=fk;
    end
end
deltas=sum(delta_ij,2)/rounds;
deltas2=sum(delta2_ij,2)/rounds;



%% Graphs %%
figure('Position', [100 100 1000 330])
hold on
yyaxis left
set(gca,'YColor', red);
plot(ns,times(:,1), '-','LineWidth', 2,'Color', red);
yyaxis right
set(gca,'YColor', blue);
plot(ns,times(:,2), '-','LineWidth', 2,'Color', blue);
legend ({'elapsed time case A', 'elapsed time case B'}) 
hold off
saveas(gcf, 'timeBsteep.png')

figure('Position', [100 100 1000 330])
hold on

yyaxis left
plot(ns,deltas, '-','LineWidth', 2,'Color', red)
plot(ns,deltas2, '-','LineWidth', 2,'Color', blue)
set(gca,'YColor', [0 0 0]);
yyaxis right
plot(ns,arrayfun(@(i) deltas(i)/ns(i),1:19), '-','LineWidth', 2,'Color', orange)
plot(ns,arrayfun(@(i) deltas2(i)/ns(i),1:19), '-','LineWidth', 2,'Color', l_blue)
legend ({'distance from minimum case A', 'distance from minimum case B', ...
    'distance from minimum/components case A', 'distance from minimum/components case B'}) 
set(gca,'YColor', [0 0 0]);
hold off
saveas(gcf, 'distBsteep.png')

%% Brown newton %%
f = @(x) fB(x);
gradf = @(x) gB(x);
Hessf = @(x) hB(x);
rho=0.5;
alpha0=1;
kmax=5000;
tolgrad=1e-8;
btmax=50;
c=1e-4;


ns=[10:10:90 100:100:1000];
len=length(ns);
rounds=1;
times=zeros(len,2);
delta_ij=zeros(len,rounds);
f_ij=zeros(len,rounds);
delta2_ij=zeros(len,rounds);
f2_ij=zeros(len,rounds);
for i=1:len
    i
    n=ns(i);
    x0 = arrayfun(@(i) 1+mod(i,2)  ,1:n);
    xmin = arrayfun(@(i) 3+log(20)/20+mod(i,2)*(3-(3+log(20)/20)),1:n);
    for j=1:rounds
        t=tolgrad*sqrt(n);
        tic
        [delta_x, fk] = Copy_of_newton_bcktrck(x0, f, gradf, Hessf, kmax, t, c, rho, btmax,xmin);
        times(i,1)=(j-1)/j*times(i,1)+1/j*toc;
        delta_ij(i,j)=delta_x;
        f_ij(i,j)=fk;
    end
    x0=zeros(n,1);
    for j=1:rounds
        t=tolgrad*sqrt(n);
        tic
        [delta_x, fk] = Copy_of_newton_bcktrck(x0, f, gradf, Hessf, kmax, t, c, rho, btmax,xmin);
        times(i,2)=(j-1)/j*times(i,2)+1/j*toc;
        delta2_ij(i,j)=delta_x;
        f2_ij(i,j)=fk;
    end
end
deltas=sum(delta_ij,2)/rounds;
deltas2=sum(delta2_ij,2)/rounds;



%% Graphs %%
figure('Position', [100 100 1000 330])
hold on
yyaxis left
set(gca,'YColor', red);
plot(ns,times(:,1), '-','LineWidth', 2,'Color', red);
yyaxis right
set(gca,'YColor', blue);
plot(ns,times(:,2), '-','LineWidth', 2,'Color', blue);
legend ({'elapsed time case A', 'elapsed time case B'}) 
hold off
saveas(gcf, 'timeBnewton.png')

figure('Position', [100 100 1000 330])
hold on

yyaxis left
plot(ns,deltas, '-','LineWidth', 2,'Color', red)
plot(ns,deltas2, '-','LineWidth', 2,'Color', blue)
set(gca,'YColor', [0 0 0]);
yyaxis right
plot(ns,arrayfun(@(i) deltas(i)/ns(i),1:19), '-','LineWidth', 2,'Color', orange)
plot(ns,arrayfun(@(i) deltas2(i)/ns(i),1:19), '-','LineWidth', 2,'Color', l_blue)
legend ({'distance from minimum case A', 'distance from minimum case B', ...
    'distance from minimum/components case A', 'distance from minimum/components case B'}) 
set(gca,'YColor', [0 0 0]);
hold off
saveas(gcf, 'distBnewton.png')

%% Functions %%
function H = hCR(x)
    n = length(x);
    princ_diag=zeros(n,1);
    princ_diag(n)=200;
    for i=1:n-1
        princ_diag(i)=-400*(-2*x(i)^2+x(i+1)-x(i)^2)+2+200;
    end
    princ_diag(1)=princ_diag(1)-200;
    up_diag=[0; -400*x(1:n-1)];
    low_diag=[-400*x(1:n-1); 0];
    H = spdiags([low_diag  princ_diag up_diag],-1:1,n,n);
end

function o = fCR(x)
    o=0;
    l=length(x);
    for i=2:l
        o=o+100*(x(i)-x(i-1)^2)^2 + (1-x(i-1))^2;
    end
end

function o = gCR(x)
    l=length(x);
    o=zeros(l,1);
    for i=2:l
        o(i-1:i)=o(i-1:i)+ ...
            [-400*(x(i)-x(i-1)^2)*x(i-1)-2*(1-x(i-1)); 200*(x(i)-x(i-1)^2)];
    end
end


function H = hB(x)
    l = length(x);
    princ_diag=zeros(l,1);
    up_diag=zeros(l,1);
    low_diag=zeros(l,1);
    k=l/2;
    for j=1:k
        i=2*j;
        princ_diag(i)=400*exp(20*(x(i-1)-x(i)));
        princ_diag(i-1)=400*exp(20*(x(i-1)-x(i)))+2+1/500;
        low_diag(i-1)=-400*exp(20*(x(i-1)-x(i)));
        up_diag(i)=-400*exp(20*(x(i-1)-x(i)));
    end
    H = spdiags([low_diag  princ_diag up_diag],-1:1,l,l);
end

function o = fB(x)
    o=0;
    t=0;
    l=length(x);
    k=l/2;
    for j=1:k
        i=2*j;
        o=o+(x(i-1)-3)^2/1000-(x(i-1)-x(i))+exp(20*(x(i-1)-x(i)));
        t=t+x(i-1)-3;
    end
    o=o+t^2;
end

function o = gB(x)
    l=length(x);
    o=zeros(l,1);
    k=l/2;
    t=0;
    for j=1:k
        i=2*j;
        t=t+x(i-1)-3;
    end
    for j=1:k
        i=2*j;
        o(i)=1-20*exp(20*(x(i-1)-x(i)));
        o(i-1)=(x(i-1)-3)/500-1+20*exp(20*(x(i-1)-x(i)))+2*t;
    end
end




