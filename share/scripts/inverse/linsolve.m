A = load('$name.gmc');
b = load('$name.imc');

x(:,1)=b(:,1);
x(:,2) = -linsolve(A,b(:,2));

save '$name_out' x '-ascii'
quit
