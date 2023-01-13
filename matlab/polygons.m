C = Polyhedron([2 2;0 2; -1 2.5]) + [2;1];
P(1) = Polyhedron([2 2;1 1; -5 3]); 
P(2) = Polyhedron([2 2;1 1; 5 -3]); 
P(3) = Polyhedron([-3 -4;1 1; -5 3]);
P(4) = Polyhedron([2 2; 6 3; 5 -3]); 

tic
N = PolyUnion(P);
Nc = N.convexHull;
Nd = (((Nc \ P(1)) \ P(2)) \ P(3)) \ P(4);
Nd = PolyUnion(Nd);
Nd = merge(Nd);
Nd = Nd.Set;

% figure(1)
% plot(Nc)
% plot(Nd,'color','green')
% plot(N,'color','yellow')
% hold on
% plot(C, 'color','cyan')

% Compute the set of translations that put C into Nc

b_i = Nc.b - max(Nc.A*C.V',[],2);
Ic = Polyhedron(Nc.A, b_i);
toc

% Compute the complement of the set of translations that put C outside Nd
for i = 1:size(Nd,1)
    b_o = Nd(i).b - min(Nd(i).A*C.V',[],2);
    O(1) = Polyhedron([Nd(i).A], b_o);

    b_o = C.b - min(C.A*Nd(i).V',[],2);
    O(2) = -Polyhedron([C.A], b_o);
    Oi(i) = intersect(O(1),O(2));
end
toc

figure(4)
plot(O)

I = Ic;
for i = 1:size(Nd,1)
    I = I \ Oi(i);
end
I = merge(PolyUnion(I));
% figure(3)
% plot(O,'color','magenta')

% hold on
% plot(Ic,'color','green')
% plot(I,'color','cyan')


% figure(2)
% plot(N,'color','yellow')
% hold on
% plot(I + C,'color','cyan')

%%
tic
Nc \ P
toc