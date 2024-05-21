close all

P = Polyhedron([7 0;9 1; 8 2]);
Q = Polyhedron([2 0; 6 6; 0 4]); 

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
hold on
plot(I,'color','green')
plot(P + Oi,'color','magenta')
plot(Q,'color','yellow')
plot(P + I,'color','orange')
plot(P, 'color', 'cyan')

save_vertices(P, 'P.txt')
save_vertices(Q, 'Q.txt')
save_vertices(I, 'I.txt')
save_vertices(Oi, 'O.txt')
save_vertices(P+I, 'P_I.txt')
save_vertices(P+Oi, 'P_O.txt')


function res = clockwise_vertices(P)
    P.minVRep;
    V = P.V;

    center = mean(V);

    angles = atan2(V(:,1)-center(1), V(:,2)-center(2));

    [~, perm] = sort(angles);

    res = V(perm,:);

end

function save_vertices(P, name)
    V = round(clockwise_vertices(P), 3);
    writematrix(V, name, 'Delimiter', ' ')
end