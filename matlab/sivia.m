f = @(x) x(1)^2 + x(2)^2 + x(1)*x(2);
Y = interval(1, 2);
X = interval([-10; -10], [10; 10]);

e = 0.1;
%%
Xc = {X};
nX = 1;


Omega = Xc;
n = nX;
% X_omeg = Polyhedron;
% for i = 1:n
%     X_omeg(i) = Omega{i}.mptPolytope.P;
% end
% X_omeg = merge(PolyUnion(X_omeg));
% figure(1)
% plot(X_omeg)
    
N = {};
nN = 0;
E = {};
nE = 0;
Xc = {};
nX = 0;

tic
while ~isempty(Omega)
    x_i = Omega{1};

    x_im = f(x_i);
    if x_im < Y
        Xc{nX+1} = x_i;
        nX = nX + 1;
    elseif isempty(intersect_interval(x_im, Y))
        N{nX+1} = x_i;
        nN = nN + 1;
    elseif max(x_i.sup - x_i.inf) < e
        E{nE+1} = x_i;
        nE = nE + 1;
    else
        [l, r] = bisect(x_i);
        Omega{n+1} = l;
        Omega{n+2} = r;
        n = n + 2;
    end
   
    
    Omega(1) = [];
    n = n - 1;
end
toc



function [l, r] = bisect(i)
    p = i.sup;
    m = i.inf;
    [~, dim] = max(p - m);
    
    offset = zeros(size(p));
    offset(dim) = (p(dim) - m(dim)) / 2;
    l = interval(m, p - offset);
    r = interval(m + offset, p);
end

function res = intersect_interval(I1, I2)
    l = max(I1.inf, I2.inf);
    h = min(I1.sup, I2.sup);
    if all(l <= h)
        res = interval(l, h);
    else
        res = [];
    end
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

