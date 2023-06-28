%%  Gram-Schmidt Orthonormalization Process

function U=GramSchmidt(V)
% Ves la matriz que contiene los vectores considerados
% almacenados en sus columnas.
[m,n] = size(V); % n vecotres, con m componentes.
U=zeros(m,n);
U(:,1)=V(:,1);
for k=2:n
    S = V(:,k);
    for j=1:k-1
        S = S-(dot(V(:,k),U(:,j))/dot(U(:,j),U(:,j)))*U(:,j);
    end
    U(:,k) = S;
end


%V1=V(1,:)/norm(V(1,:))