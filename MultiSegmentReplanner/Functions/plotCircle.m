function plotCircle(xc, yc, r)
    th = 0:pi/50:2*pi;
    x = r * cos(th) + xc;
    y = r * sin(th) + yc;
    plot(x, y);
end