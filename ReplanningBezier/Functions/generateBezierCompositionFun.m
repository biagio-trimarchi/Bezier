function generateBezierCompositionFun(m, n)
%% Symbolic Composition
    % Curve Orders

    P = sym('p', [m+1, 1]);
    U = sym('u', [n+1, 1]);

    A = P;

    for kk = 1:m
        for ii = 0:m-kk
            for jj = 0:kk*n
                iIdx = ii+1;
                kIdx = kk+1;
                jIdx = jj+1;

                sum = 0;
                for ll = max([0, jj-n]):min([jj, kk*n-n])
                    lIdx = ll+1;

                    w1 = binomial(kk*n-n, ll);
                    w2 = binomial(n, jj-ll);
                    w3 = ((1-U(jIdx-lIdx+1,1))*A(iIdx, kIdx-1, lIdx)) + ...
                     (U(jIdx-lIdx+1,1)*A(iIdx+1, kIdx-1, lIdx));

                    sum = sum + (w1*w2*w3);
                end

                A(iIdx, kIdx, jIdx) = (1/binomial(kk*n, jj))*sum;
            end
        end
    end

    C(:,1) = A(1, m+1, 1:((m*n)+1));
    C = simplify(C); % New Control Points
    matlabFunction(C,'File','./Functions/BezierComposition');
end