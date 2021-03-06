function plotWall()
    alphaFace = 0.4;
    n = 16;
    wall(1).coord = [5, 5, 8, 8;
                     0, 0, 0, 0;
                     0, 9, 9, 0];
      
    wall(2).coord = [5, 5, 5, 5;
                     0, 0, 3, 3;
                     0, 9, 9, 0];

    wall(3).coord = [5, 5, 8, 8;
                     3, 3, 3, 3;
                    3, 6, 6, 3];

    wall(4).coord = [8, 8, 8, 8;
                     0, 0, 3, 3;
                     0, 9, 9, 0];
                
    wall(5).coord = [5, 5, 5, 5;
                     3, 3, 6, 6;
                     0, 3, 3, 0];
             
    wall(6).coord = [8, 8, 8, 8;
                     3, 3, 6, 6;
                     0, 3, 3, 0];

    wall(7).coord = [5, 5, 5, 5;
                     3, 3, 6, 6;
                     6, 9, 9, 6];
             
    wall(8).coord = [8, 8, 8, 8;
                     3, 3, 6, 6;
                     6, 9, 9, 6];

    wall(9).coord = [5, 5, 8, 8;
                     6, 6, 6, 6;
                     3, 6, 6, 3];
      
    wall(10).coord = [5, 5, 5, 5;
                      6, 6, 9, 9;
                      0, 9, 9, 0];

    wall(11).coord = [5, 5, 8, 8;
                      9, 9, 9, 9;
                      0, 9, 9, 0];

    wall(12).coord = [8, 8, 8, 8;
                      6, 6, 9, 9;
                      0, 9, 9, 0];
              
    wall(13).coord = [5, 5, 8, 8;
                      3, 6, 6, 3;
                      3, 3, 3, 3];

    wall(14).coord = [5, 5, 8, 8;
                      3, 6, 6, 3;
                      6, 6, 6, 6];

    wall(15).coord = [5, 5, 8, 8;
                      0, 9, 9, 0;
                      0, 0, 0, 0];
              
    wall(16).coord = [5, 5, 8, 8;
                      0, 9, 9, 0;
                      9, 9, 9, 9];

    for i =1:n
        patch(wall(i).coord(1, :), wall(i).coord(2, :), wall(i).coord(3, :), 'black', 'FaceAlpha', alphaFace, 'EdgeColor', 'none')
    end
end