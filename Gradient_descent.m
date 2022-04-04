%x0 = [15 15]'; % Zadavanje početne točke
x0 = [-9 -9]';% Za Rosenbrockovu funkciju
f = @(a,b) -(10-a).^2+50*(10*b-a.^2).^2 ;
%f = @(x1,x2) x1.^2 + x1.*x2 + 3*x2.^2; % Tipična kvadratna funkcija
% redefine objective function syntax for use with optimization:
f2 = @(x) f(x(1),x(2));

err = 0.5e-9; 
errg = 0.5e-3;
maxiter=1000;
xnorm = inf;
gnorm = inf;
x = x0;
niter = 0;
figure(1);
a = -10:0.05:10; % Za Rosenbrockovu funkciju 
b = -10:0.05:10; % Za Rosenbrockovu funkciju
%a = -20:0.05:20; 
%b = -20:0.05:20;
[X,Y] = meshgrid(a,b);
Z = -(10-X).^2+50*((10*Y-(X.^2)).^2) ;
%Z = X.^2+X.*Y+3*Y.^2;
surf(X,Y,Z)
colormap(jet(31));
shading interp;
hold on

alfa = input('Za unaprijd zadani alfa unesi 1,\nza računanje alfa pomoću hesijana 2,\na za računanje alfa pomoću procjenjenih vrijednosti 3');

switch alfa
    case 0
        %alfa=0.01;
        alfa=0.00001; % Za Rosenbrockovu funkciju
        flag = 0;
    case 1
        flag = 1;
    case 2
        fold= f2(x0);
        %alfa=0.3;
        alfa=3e-5; % Za Rosenbrockovu funkciju
        flag = 2;
    otherwise
        disp('Niste unijeli 0,1 ili 2.')
        error('Alfa nije definiran.')
end

while and(gnorm >= errg, and(niter <= maxiter, xnorm >= err))
    gnormp=gnorm;
    g = grad(x);
    if (flag == 1)
        h = hessian(x);
        alfa = (g'*g)/(g'*h*g);
    end
  
    % take step:
    xnew = x - alfa*g;
    if (flag == 2)
        fnew = f2(xnew);
        alfa = -((g)'*g*alfa^2)/(2*(fnew-fold-alfa*(g)'*g));
        % if alfa < 1e-7  % Nužno kod Rosenbrockove funkcije za
%         približavanje minimumu
        %     alfa = 1e-7;
        % end
        fold = fnew;
  
%          if (gnorm-gnormp)<0.01 %Ubrzava za kvadratnu funkciju
%             flag=3;
%             alfa=0.25;
%             err = 0.5e-3;
%          end
        
    end
    xnorm = norm(x-xnew);
    gnorm = norm(g);
    % update termination metrics
    niter = niter + 1;
    
    % check step
    if ~isfinite(xnew)
        display(['Number of iterations: ' num2str(niter)])
        error('x is inf or NaN')
    end
    % plot current point
    if (flag == 0)
    plot3([x(1) xnew(1)],[x(2) xnew(2)],[f(x(1),x(2)) f(xnew(1),xnew(2))],'ro-')
    refresh
    end
    if(flag == 1)
    plot3([x(1) xnew(1)],[x(2) xnew(2)],[f(x(1),x(2)) f(xnew(1),xnew(2))],'g*-')  
    end
    if or(flag == 2 , flag == 3)
    plot3([x(1) xnew(1)],[x(2) xnew(2)],[f(x(1),x(2)) f(xnew(1),xnew(2))],'yx-')
    end
  x = xnew;
end
xopt = x;
fopt = f2(xopt);
niter = niter - 1;
fprintf('Minimum je nadjen na poziciji x = %.4f, y = %.4f',xopt(1),xopt(2))
fprintf('\nVrijednost funkcije u minimumu je %.4f', fopt)
fprintf('\nBroj iteracija %d',niter)
function g = grad(x)
g = [2*(10-x(1))-200*(10*x(2)*x(1)-x(1)^3) ; 1000*(10*x(2)-x(1)^2)];
%g = [2*x(1) + x(2);x(1) + 6*x(2)];
end
function h = hessian(x)
h=[600*x(1)^2-2000*x(2)-2 -2000*x(1);-2000*x(1) 10000];
%h=[2, 1;1 6];
end
