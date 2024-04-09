% Natural Frequency (OMEGA) and Mode Shape (PHI) of a system
% [M] and [K] should be defined by user in the next section.
% --------------------------------------------
M=[1 0; 0 1];
K=[2 -1;-1 1];
% --------------------------------------------
[V,D]=eig(K,M);
[W,KK]=sort(diag(D));
V=V(:,KK);
PHI=V*inv(sqrt(diag(diag(V'*M*V))))  % Mode Shape
OMEGA=diag(sqrt(PHI'*K*PHI))         % Natural Frequency
% --------------------------------------------