function A = Jacobian_computation_parametric(x1, y1, x6, y6, w1, w6, w9, w15, X_0,points)
    syms x1 y1 x6 y6 x9 y9 x15 y15 w1 w6 w9 w15
    
    x9=points(3,3);
    y9=points(3,2);
    x15=points(4,3);
    y15=points(4,2);
    
    % Distances
    L_sym(1) = sqrt((x1-x6)^2+(y1-y6)^2);
    L_sym(2) = sqrt((x1-x9)^2+(y1-y9)^2);
    L_sym(3) = sqrt((x6-x9)^2+(y6-y9)^2);
    L_sym(4) = sqrt((x1-x15)^2+(y1-y15)^2);
    L_sym(5) = sqrt((x9-x15)^2+(y9-y15)^2);
    % Directions
    L_sym(6) = atan((y6-y1)/(x6-x1))-w1;%check quadron
    L_sym(7) = atan((y15-y1)/(x15-x1))-w1;
    L_sym(8) = atan((y1-y6)/(x1-x6))-w6;
    L_sym(9) = atan((y9-y6)/(x9-x6))-w6; 
    L_sym(10) = atan((y15-y9)/(x15-x9))-w9;
    L_sym(11) = atan((y1-y9)/(x1-x9))-w9; 
    L_sym(12) = atan((y6-y9)/(x6-x9))-w9;
    L_sym(13) = atan((y1-y15)/(x1-x15))-w15;
    L_sym(14) = atan((y9-y15)/(x9-x15))-w15; 

    J = jacobian(L_sym,[x1, y1, x6, y6, w1, w6, w9, w15]);
    A = subs(J, [x1, y1, x6, y6, w1, w6, w9, w15],X_0');
end
