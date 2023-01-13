function Xbar = cz_prediction(func,fp,X0,u,flagH)

if size(X0.Z,2)-1 > 3
    ng_desired = 3;
    nc_desired = 2;
    n_order    = size(X0.Z,1);
    order = (ng_desired-nc_desired)/n_order;
    
    Xt = reduce(X0,'girard',order,nc_desired);
    
    if ~isempty(Xt)
        X0 = Xt;
    else
        nc_desired = 1;
        order = (ng_desired-nc_desired)/n_order;
        X0 = reduce(X0,'girard',order,nc_desired);
    end
end

% X0 = conZonotope(removeRedundancies(mptPolytope(X0)));
% X0 = reduce(X0,'redConstr');

c = X0.Z(:,1);
G = X0.Z(:,2:end);
A = X0.A;
b = X0.b;

%% f(eta) and its Bounds

eta = sym('eta',[size(G,2) 1]);

if isempty(X0.A)
    inv_eta = interval(-ones(size(G,2),1),ones(size(G,2),1));
else
    ss=size(size(pinv(A),1),1);
    rs=zeros(ss,1);
    II=eye(ss);
    AAf=pinv(A)*A;
    MM=Inf;
    for i=1:ss
        if II(:,i)~=AAf(:,i)
            rs(i)=1;
        end
    end
    inv_eta = interval(max(-ones(size(G,2),1),pinv(A)*b-MM*rs),min(ones(size(G,2),1),pinv(A)*b+MM*rs));
end

%%

func_xi = @(eta) func(c+G*eta,u,fp);
func_xi = func_xi(eta);

func_xis = cell(1,size(func_xi,1));

for m = 1:size(func_xi,1)
    func_xis{m} = matlabFunction(func_xi(m,:),'Vars',{eta});
end

%%

[HH,fupp,floww]  = decomposition_bounds(func_xis,eta,inv_eta);
[H,gd_up,gd_low,~,~] = Hcombinations(HH,fupp,floww,'P',flagH);

%% Constrained Zonotopes

cf = cell(1,length(H));
Gf = cell(1,length(H));
Af = cell(1,length(H));
bf = cell(1,length(H));
cz = cell(1,length(H));

for k = 1%:length(H)
    
    cf{k} = 1/2*(gd_low{k} + gd_up{k});
    Gf{k} = [H{k} 1/2*diag(gd_up{k} - gd_low{k})];
    Af{k} = [A zeros(size(A,1),size(Gf{k},2)-size(A,2))];
    bf{k} = b;
    
    cz{k} = conZonotope(cf{k},Gf{k},Af{k},bf{k});
    
end

Xbar = cz{1};

return;

czi = mptPolytope(cz{1});
for i = 2:length(cz)
    czi = removeRedundancies(and(mptPolytope(cz{i}),czi));
end

Xbar = conZonotope(czi);