function [J, dJ] = cost2(x, n, d, Cd)
   J = (x-Cd)*(x-Cd)';
   dJ = (x-Cd);
end