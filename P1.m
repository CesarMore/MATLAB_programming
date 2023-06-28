%% Tarea: calcular las matrices expotenciales
syms a1 a2 s
A = [0,1;-a2,-a1]
X = inv(s*eye(2)-A)
%% caso 1
A1 = subs(A,[a1 a2], [2 1])
X1s = subs(X,[a1 a2], [2 1])
X1t = ilaplace(X1s)
%% caso 2
A2 = subs(A,[a1 a2], [0 0])
X2s = subs(X,[a1 a2], [0 0])
X2t = ilaplace(X2s)
%% caso 3
A3 = subs(A,[a1 a2], [0 1])
X3s = subs(X,[a1 a2], [0 1])
X3t = ilaplace(X3s)
%%
[A1, A2, A3]
[X1t, X2t, X3t]
