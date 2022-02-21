function plotSphere(xc, yc, zc, r)
    [X, Y, Z] = sphere();
    X = X*r;
    Y = Y*r;
    Z = Z*r;
    surf(X+xc, Y+yc, Z+zc, 'FaceAlpha', 0.1, 'EdgeColor', 'none')
end