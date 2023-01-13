x_c = 70
r_c = 0.005
x_l = 3
r_l = 0.05
r_0 = 1
v_s = 1

A_1 = [-r_l/x_l 0; 0 -1/(x_c*(r_c+r_0))]
A_2 = [-1/x_l*(r_l + r_0*r_c/(r_c+r_0)) -r_0/(x_l*(r_c+r_0)); r_0/(x_c*(r_c+r_0)) -1/(x_c*(r_c+r_0))]
b = [v_s/x_l; 0]

t_s = 0.5

sys1 = ss(A_1, b, zeros(1,2),0)
sys1d = c2d(sys1, 0.5, 'zoh')

sys2 = ss(A_2, b, zeros(1,2),0)
sys2d = c2d(sys2, 0.5, 'zoh')