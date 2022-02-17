function [x,pleft,pright]=castel(p,t)
% The de Casteljau algorithm.
% Input is the control points given in a row vector p and the desired time t.
% Output is the value of the Bezier curve at time t, and also the points
% giving the left and right part of the subdivision for the given t, given in row vectors.
n=size(p,2); %n is the number of points
dim=size(p,1); %the dimension of the points, i.e. 2d or 3d
p=p'; %transpose the matrix p. Now, the 1st column is the x coordinate, second is the y, 3rd is the z.

%initiate control polygons
pleft=[p(1,:)'];
pright=[p(end,:)'];
if dim==2 % check if the points are 2d or 3d
    for i=1:n-1 % careful, here n refers to the number of points. As in most books, for an nth degree curve, there are n+1 points p 0,...,p n.
        for j=1:n-i % again here we don't start counting from zero but from 1, up to n-i
            p(j,1)=(1-t)*p(j,1)+t*p(j+1,1); % x coordinate
            p(j,2)=(1-t)*p(j,2)+t*p(j+1,2); % y coordinate
        end
        
        % I also want to keep the control polygons giving the subdivision
        pleft=[pleft, p(1,:)']; %each new point is placed last in the list
        pright=[p(end-i,:)',pright]; %each new point is placed first in the list
    end
    x=[p(1,1);p(1,2)];
else
    for i=1:n-1
        for j=1:n-i
            p(j,1)=(1-t)*p(j,1)+t*p(j+1,1); % x coordinate
            p(j,2)=(1-t)*p(j,2)+t*p(j+1,2); % y coordinate
            p(j,3)=(1-t)*p(j,3)+t*p(j+1,3); % z coordinate
        end
        % I also want to keep the control polygons giving the subdivision
        pleft=[pleft, p(1,:)']; %each new point is placed last in the list
        pright=[p(end-i,:)',pright]; %each new point is placed first in the list
    end
    x=[p(1,1);p(1,2);p(1,3)];
end