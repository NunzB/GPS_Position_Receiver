function [test,coord_graus] = isValidcoordinates(coordinates)

if (length(coordinates)~=10)
    test=0;
    coord_graus=0;
    return;
end

if (coordinates(3)~='-'||coordinates(6)~='-'||coordinates(9)~='-')
    test=0;
    coord_graus=0;
    return;
end

dd = str2double(coordinates(1:2));
mm = str2double(coordinates(4:5));
ss = str2double(coordinates(7:8));
N_S= coordinates(10);

if N_S == 'N' || N_S =='E'
    coord_graus = DMStoD ([dd mm ss]);
    test=1;
    return
else if N_S == 'S' || N_S =='W'
    coord_graus = -DMStoD ([dd mm ss]);   
    test=1;
    return
    end
end

if (abs(coord_graus)>180)
if (abs(coord_graus)>90)
    test=0;
    coord_graus=0;
    return;
end
end
end

