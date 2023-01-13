function [A,B,C,D,Jphi_x_upp,Jphi_w_upp,Jpsi_x_upp,Jpsi_v_upp] = FYM_gains_noisy_3_compare(f,h,interval_X,interval_W,interval_V)

solver=sdpsettings('solver', 'sedumi', 'sedumi.eps', 1e-8, ...
                'sedumi.cg.qprec', 1, 'sedumi.cg.maxiter', 49, ...
                'sedumi.stepdif', 2);
nx=size(interval_X,1);
nw=size(interval_W,1);
nv=size(interval_V,1);
%if interval_V~=[]
if isempty(interval_V)==0
p = size(h(interval_X,interval_V),1);
else
    p = size(h(interval_X),1);
end

x = sym('x',[nx 1]);
w = sym('w',[nw 1]);
v = sym('v',[nv 1]);

grad_f_x = jacobian(f(x,w),x);
grad_f_w = jacobian(f(x,w),w);
%if interval_V~=[]
if isempty(interval_V)==0
grad_h_x = jacobian(h(x,v),x);
grad_h_v = jacobian(h(x,v),v);
else
grad_h_x = jacobian(h(x),x);
%grad_h_v = jacobian(h(x,v),v);
end

gradfunc_f_x = matlabFunction(grad_f_x,'Vars',{x,w});   
gradfunc_f_w = matlabFunction(grad_f_w,'Vars',{x,w}); 
%if interval_V~=[]
if isempty(interval_V)==0
gradfunc_h_x = matlabFunction(grad_h_x,'Vars',{x,v});   
gradfunc_h_v = matlabFunction(grad_h_v,'Vars',{x,v}); 
else
gradfunc_h_x = matlabFunction(grad_h_x,'Vars',{x});   
%gradfunc_h_v = matlabFunction(grad_h_v,'Vars',{x,v});
end

Jfx= gradfunc_f_x(interval_X,interval_W);
Jfw= gradfunc_f_w(interval_X,interval_W);
%if interval_V~=[]
if isempty(interval_V)==0
Jhx= gradfunc_h_x(interval_X,interval_V);
Jhv= gradfunc_h_v(interval_X,interval_V);
else
    Jhx= gradfunc_h_x(interval_X);
end

Jfx_upp = interval(Jfx).sup;
Jfx_low = interval(Jfx).inf;

Jfw_upp = interval(Jfw).sup;
Jfw_low = interval(Jfw).inf;

Jhx_upp = interval(Jhx).sup;
Jhx_low = interval(Jhx).inf;

%if interval_V~=[]
if isempty(interval_V)==0
Jhv_upp = interval(Jhv).sup;
Jhv_low = interval(Jhv).inf;
end

%Jy_upp = [1 0];
%Jy_low = [1 0];

        A=zeros(nx,nx);
        C=zeros(p,nx);
        B=zeros(nx,nw);
        D=zeros(p,nv);
        for i=1:nx
            for j=1:nx
                if Jfx_low(i,j)>=0 || Jfx_upp(i,j)<=0
                    A(i,j)=0;
                elseif abs(Jfx_upp(i,j))>abs(Jfx_low(i,j))
                    A(i,j)=Jfx_low(i,j);
                  %  H(i,j)=Jf_upp(i,j);% test
                else
                    A(i,j)=Jfx_upp(i,j);
                  %  H(i,j)=Jf_low(i,j);%test
                end
                % test
               % if A(i,j)~=0
               % A(i,j)=Jfx_upp(i,j);
               % end
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            end
                
            for k=1:nw 
                if Jfw_low(i,k)>=0 || Jfw_upp(i,k)<=0
                    B(i,k)=0;
                elseif abs(Jfw_upp(i,k))>abs(Jfw_low(i,k))
                    B(i,k)=Jfw_low(i,k);
                  %  H(i,j)=Jf_upp(i,j);% test
                else
                    B(i,k)=Jfw_upp(i,k);
                  %  H(i,j)=Jf_low(i,j);%test
                end
                % test
              %  if B(i,k)~=0
              %  B(i,k)=Jfw_upp(i,k);
              %  end
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            end
        end
        
        for i=1:p
            for j=1:nx
                if Jhx_low(i,j)>=0 || Jhx_upp(i,j)<=0
                    C(i,j)=0;
                elseif abs(Jhx_upp(i,j))>abs(Jhx_low(i,j))
                    C(i,j)=Jhx_low(i,j);
                  %  H(i,j)=Jf_upp(i,j);% test
                else
                    C(i,j)=Jhx_upp(i,j);
                  %  H(i,j)=Jf_low(i,j);%test
                end
                % test
              %  if C(i,j)~=0
              %  C(i,j)=Jhx_upp(i,j);
              %  end
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            end
     %if interval_V~=[] 
     if isempty(interval_V)==0
            for k=1:nv 
                if Jhv_low(i,k)>=0 || Jhv_upp(i,k)<=0
                    D(i,k)=0;
                elseif abs(Jhv_upp(i,k))>abs(Jhv_low(i,k))
                    D(i,k)=Jhv_low(i,k);
                  %  H(i,j)=Jf_upp(i,j);% test
                else
                    D(i,k)=Jhv_upp(i,k);
                  %  H(i,j)=Jf_low(i,j);%test
                end
                % test
              %  if D(i,k)~=0
              %  D(i,k)=Jhv_upp(i,k);
              %  end
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            end
     end
        end
  %A=H_f;C  
%test:
%H=Jf_upp;
%H=Jf_low;
%%%%%%%%%%%%%%%%%%%%%%%%
J_phi_x=Jfx-A;
J_phi_w=Jfw-B;
J_psi_x=Jhx-C;
%J_psi_v=Jhv-D;
Jphi_x_upp = interval(J_phi_x).sup;
Jphi_w_upp = interval(J_phi_w).sup;
Jphi_x_low = interval(J_phi_x).inf;
Jphi_w_low = interval(J_phi_w).inf;
Jpsi_x_upp = interval(J_psi_x).sup;
Jpsi_x_low = interval(J_psi_x).inf;
%Jpsi_v_upp = interval(J_psi_v).sup;
%Jpsi_v_low = interval(J_psi_v).inf;
%if interval_V~=[]
if isempty(interval_V)==0
J_psi_v=Jhv-D;
Jpsi_v_upp = interval(J_psi_v).sup;
Jpsi_v_low = interval(J_psi_v).inf;
else
    J_psi_v=[];
    Jpsi_v_upp = [];
    Jpsi_v_low = [];
end
%Jphi_upp = J_phi.sup;
%Jphi_low = J_phi.inf;

%%
%F_upp_phi = 2*max(Jphi_upp, zeros(n,n)) - Jphi_low ;
%F_upp_ksi = 2*max(Jy_upp-C, zeros(p,n)) - Jy_low + C;

% test:
F_upp_phi_x = Jfx_upp- Jfx_low ;
F_upp_phi_w = Jfw_upp- Jfw_low ;
F_upp_psi_x = Jhx_upp- Jhx_low ;
%if interval_V~=[]
if isempty(interval_V)==0
F_upp_psi_v = Jhv_upp- Jhv_low ;
end

F_upp_phi_x = 2*max(Jphi_x_upp, zeros(nx,nx)) - Jphi_x_low ;
F_upp_phi_w = 2*max(Jphi_w_upp, zeros(nw,nw)) - Jphi_w_low ;
%F_upp_ksi = 2*max(Jy_upp-C, zeros(p,n)) - Jy_low + C;

F_upp_phi_x = max(Jphi_x_upp, zeros(nx,nx)) - min(Jphi_x_low, zeros(nx,nx)) ;
F_upp_phi_w = max(Jphi_w_upp, zeros(nw,nw)) - min(Jphi_w_low, zeros(nw,nw)) ;
F_upp_psi_x = max(Jpsi_x_upp, zeros(size(Jpsi_x_upp))) - min(Jpsi_x_low, zeros(size(Jpsi_x_low))) ;
if isempty(interval_V)==0
F_upp_psi_v = max(Jpsi_v_upp, zeros(size(Jpsi_v_upp))) - min(Jpsi_v_low, zeros(size(Jpsi_v_low))) ;
end

% test:
%F_upp_phi_x = Jfx_upp- Jfx_low ;
%F_upp_phi_w = Jfw_upp- Jfw_low ;
%F_upp_ksi_x = Jhx_upp- Jhx_low ;
%if interval_V~=[]
%F_upp_ksi_v = Jhv_upp- Jhv_low ;
%end

