function C = fastintersect(A,B)
P = zeros(1, max(max(A),max(B)) ) ;
P(A) = 1;
C = B(logical(P(B)));
