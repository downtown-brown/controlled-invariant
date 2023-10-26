close all

P = Polyhedron([0 0;2 1; 1 2]);
Q = Polyhedron([4 0; 8 6; 2 4]); 

figure(1)
plot(P, 'color', 'cyan')
hold on
plot(Q, 'color', 'yellow')

%%
b_i = Q.b - max(Q.A*P.V',[],2);
I = Polyhedron(Q.A, b_i);


% Compute the complement of the set of translations that put C outside Nd
b_o = Q.b - min(Q.A*P.V',[],2);
O(1) = Polyhedron([Q.A], b_o);

b_o = P.b - min(P.A*Q.V',[],2);
O(2) = -Polyhedron([P.A], b_o);
Oi = intersect(O(1),O(2));

figure(4)
plot(Oi)

figure(3)
plot(I,'color','green')

figure(2)
plot(Q,'color','yellow')
hold on
plot(P + I,'color','cyan')

figure(12)
plot(P + Oi,'color','cyan')
hold on
plot(Q,'color','yellow')
