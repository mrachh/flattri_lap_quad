
function A = Afun_lap_te(i,j,x,nu,area,P,S,contrast,wuse)
% AFUN(I,J) computes entries of the matrix A to be factorized at the
% index sets I and J.  This handles the near-field correction.
if isempty(i) || isempty(j)
  A = zeros(length(i),length(j));
  return
end
[I,J] = ndgrid(i,j);
A = 2*bsxfun(@times,lap_neumann_kernel(x(:,i),x(:,j),nu(:,i)),area(j));
A = bsxfun(@times,contrast(i),A);
M = -2*spget_quadcorr(i,j,P,S);
idx = abs(M) ~= 0;
A(idx) = A(idx) + M(idx);
A(I == J) = A(I == J) + 1.0;
A = A + wuse*ones(length(i),1)*area(j);

A = bsxfun(@times,sqrt(area(i)).',A);
A = bsxfun(@times,A,1.0./sqrt(area(j)));
end
