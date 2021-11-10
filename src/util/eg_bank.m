function [X,V,name] = eg_bank(sides)
%list of embedding parameters for standard quasi-regular polygons
    switch sides
        case 1
            error('No one sided shape. If you meant the screen, choose n=2');
        case 2 %screen
            V =  [0 0; 1 0];
            X = struct('N',2,'M',2,'p',1,'q',2,'V',V);
            name = 'screen';
        case 3 %triangle
            V = [0.5 sqrt(3)/2; 1 0; 0 0];
            X = struct('N',3,'M',12,'p',3,'q',5,'V',V);
            name = 'equilateral triangle';
        case 4 %square
            V = sqrt(2)*[-.5 -.5; -.5 .5; .5 .5; .5 -.5];
            X = struct('N',4,'M',8,'p',2,'q',3,'V',V);
            name = 'square';
        case 5 %pentagon
            % work out the vertices later
            X = struct('N',5,'M',30,'p',5,'q',7);
            name = 'pentagon';
        case 6 % hexagon
            V = [3/2 sqrt(3)/2; 1 0; 0 0; -1/2 sqrt(3)/2; 0 sqrt(3); 1 sqrt(3)];
            X = struct('N',6,'M',18,'p',3,'q',4,'V',V);
            name = 'hexagon';
        otherwise
            error('Only go as high as n=6 for now');
    end
    X.twoPhi = X.q*pi/X.p;
end

