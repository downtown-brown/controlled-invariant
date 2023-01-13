C = Polyhedron([2 2;0 2; -1 2.5]) + [2;1];

N = Polyhedron([-1 1 1; 1 1 1; 1 -1 1; -1 -1 1; -1 -1 -1; 1 -1 -1; 1 1 -1; -1 1 -1])*10;
T = [1 0 0.1;0 1 0.1];
t = [0; 1];

C_t = myInvAffineMap(C, T, t);

C_new = T*C_t + t;

W_in = randn(10, 2);
b_in = randn(10, 1);
W_1 = randn(10,10);
b_1 = randn(10,1);
W_out = randn(1, 10);
b_out = randn(1, 1);

f = @(x, y) W_out*(max((W_1*(W_in*[x;y] + b_in) + b_1), 0)) + b_out;

N2 = Polyhedron(A = [eye(2); -eye(2)], b = [1; 1; 1; 1])
M = Polyhedron(A=-eye(10), b = zeros(10,1))

L1 = Polyhedron(A=N2.A*W_in', b = N2.b*W_in' + b_in);
% L1 = W_in*N2 + b_in

L2 = W_1*L1 + b_1
L3 = intersect(L2, M)

L4 = W_out*L3 + b_out

[X, Y] = meshgrid(linspace(-1,1,100), linspace(-1,1,100));
Z = zeros(100,100);
for i = 1:100
    for j = 1:100
        Z(i,j) = f(X(i,j), Y(i,j));
    end
end

surf(X, Y, Z)

%%
P = Polyhedron([1 -.1; -10 3; 0 2; 2 0]);
N2 = Polyhedron(A = [eye(2); -eye(2)], b = [1; 1; 1; 1])*10;

Pre = [];

for v1 = [0,1]
    for v2 = [0,1]
        v = [v1; v2];
        tmp = Polyhedron([P.A*diag(v); -diag(v); diag(1-v)], [P.b; zeros(2,1); zeros(2,1)]);
        tmp = intersect(tmp, N2);
        Pre = [Pre tmp];
    end
end

plot(Pre)

