function RIF = RIF(Q,P,A)
%RIF = RIF(Q,P,A)
%
%RIF returns the recursive impact factor for a set of journals
%P is a KxJ matrix. The rows have the number of citations by paper k to
%papers published in journal j.
%Q is a JxK matrix, mapping paper k to the journal j in which it was
%published.
%A is a diagonal JxJ matrix with the number of papers published in journal
%j. (I know we could construct A from Q.)

C = Q*P;
D = diag(sum(C,1));
M = inv(A)*C*inv(D)*A;
[V,D] = eig(M);

RIF = V(:,1);
end