function temp = temperature2d(i, j, T, BC, w, e, s, n, source_u)
numerator = 0;
switch BC
    case 'nw'
        numerator = e*T(i,j+1) + s*T(i+1,j) + source_u;
    case 'sw'
        numerator = e*T(i,j+1) + n*T(i-1,j) + source_u;
    case 'ne'
        numerator = w*T(i,j-1) + s*T(i+1,j) + source_u;
    case 'se'
        numerator = w*T(i,j-1) + n*T(i-1,j) + source_u;
    case 'w'    
        numerator = e*T(i,j+1) + s*T(i+1,j) + n*T(i-1,j) + source_u;
    case 'e'
        numerator = w*T(i,j-1) + s*T(i+1,j) + n*T(i-1,j) + source_u;
    case 'n'
        numerator = e*T(i,j+1) + w*T(i,j-1) + s*T(i+1,j) + source_u;
    case 's'
        numerator = e*T(i,j+1) + w*T(i,j-1) + n*T(i-1,j) + source_u;
    case 'interior'
        numerator = e*T(i,j+1) + w*T(i,j-1) + s*T(i+1,j) + n*T(i-1,j) + source_u;
end
denominator = e + w + s + n;
temp = numerator/denominator;
end