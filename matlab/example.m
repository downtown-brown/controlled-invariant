% Example from
% https://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=7084969

% n = 2;
% mu = 0.5;
% dt = 0.01;
% 
% u_i = interval(-2,2);
% 
% f1 = @(x) x + dt.*[x(2); x(1)];
%         
% g = @(x) dt*[(1-mu)*x(1); 4*(1-mu)*x(2)];
% 
% f = @(x,u) f1(x) + g(x).*u;
% 
% X = interval([-1;-1],[1;1]);
% U = interval(-2,2);
% % for i = 1:N
% %     while running
% %         
% %     end
% % end
% 
% % [A,B,C,D,Jphi_x_upp,Jphi_w_upp,Jpsi_x_upp,Jpsi_v_upp] = FYM_gains_noisy_3_compare(f,@(x,v) 0,X,0,U)
% A = [1 dt; dt 1];
% B = @(x) dt*[(1-mu)*x(1); 4*(1-mu)*x(2)];
% phi = @(X,U,x) U*dt*[(1-mu)*(X(1)-x(1)); 4*(1-mu)*(X(2)-x(2))];
% 
% e = .001;

n = 2;
dt = 0.01;

m = 0.2;
g = 9.8;
l = 0.3;
J = 0.006;
b = 0.1;

A = [1 dt; m*g*l/J*dt (1 - dt*b/J)];
B = @(x) dt*[0; l/J*cos(x(1))];

phi = @(X, U, x) dt*U*[0; l/J*(cos(X(1)) - cos(x(1)))];

e = 0.001;

X = interval([-0.05;-0.01], [0.05;0.01]);   
U = interval(-0.1,0.1);


%%
Xc = {X};
nX = 1;

for t = 1:100

Omega = Xc;
n = nX;
X_omeg = Polyhedron;
for i = 1:n
    X_omeg(i) = Omega{i}.mptPolytope.P;
end
X_omeg = merge(PolyUnion(X_omeg));
figure(1)
plot(X_omeg)
    
N = {};
nN = 0;
E = {};
nE = 0;
Xc = {};
nX = 0;
i = 0;
tic
while ~isempty(Omega)
    x_i = Omega{1};
    %fprintf("Considering interval [%f, %f] x [%f, %f]\n", x_i.inf(1), x_i.sup(1), x_i.inf(2), x_i.sup(2)) 
    
    x_m = (x_i.sup + x_i.inf)/2;
    tmp = phi(x_i,U,x_m);
    Xo = A*x_i.mptPolytope.P + B(x_m)*U.mptPolytope.P + tmp.mptPolytope.P;
    
    Xz = A*x_i.mptPolytope.P + tmp.mptPolytope.P;
    
    Xi = intersect_polyunion(X_omeg, Xo);
    if isempty(Xi.Set)
        %fprintf("rejected\n")
        N{nN+1} = x_i;
        nN = nN + 1;
    elseif check_overlap(Xi, Xz)
        %fprintf("accepted\n")
        Xc{nX+1} = x_i;
        nX = nX + 1;
    elseif max(x_i.sup - x_i.inf) < e
        %fprintf("boundary\n")
        E{nE+1} = x_i;
        nE = nE + 1;
    else
        %fprintf("bisected\n")
        [l, r] = bisect(x_i);
        Omega{n+1} = l;
        Omega{n+2} = r;
        n = n + 2;
    end
   
    
    Omega(1) = [];
    n = n - 1;
    i = i + 1;
end
fprintf("Considered %d intervals\n", i);
toc
end


function [l, r] = bisect(i)
    p = i.sup;
    m = i.inf;
    [~, dim] = max(p - m);
    
    offset = zeros(size(p));
    offset(dim) = (p(dim) - m(dim)) / 2;
    l = interval(m, p - offset);
    r = interval(m + offset, p);
end

function res = intersect_polyunion(Xu, Xp)
    res = Polyhedron;
    for i = 1:Xu.Num
        res(i) = intersect(Xu.Set(i), Xp);
    end
    res = merge(PolyUnion(res));
        
end

function res = check_overlap(N, C)
% Check if C fits inside N
Nc = N.convexHull;

for i = 1:N.Num
    Nd = Nc \ N.Set(i);
end
Nd = PolyUnion(Nd);
Nd = merge(Nd);
Nd = Nd.Set;

% Compute the set of translations that put C into Nc

b_i = Nc.b - max(Nc.A*C.V',[],2);
Ic = Polyhedron(Nc.A, b_i);

% Compute the complement of the set of translations that put C outside Nd
for i = 1:size(Nd,1)
    b_o = Nd(i).b - min(Nd(i).A*C.V',[],2);
    O(1) = Polyhedron([Nd(i).A], b_o);

    b_o = C.b - min(C.A*Nd(i).V',[],2);
    O(2) = -Polyhedron([C.A], b_o);
    Oi(i) = intersect(O(1),O(2));
end

I = Ic;
for i = 1:size(Nd,1)
    I = I \ Oi(i);
end
I = merge(PolyUnion(I));
res = ~isempty(I.Set);

end

