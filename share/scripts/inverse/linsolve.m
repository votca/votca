A = load('%NAME.gmc');
b = load('%NAME_m.imc');

b(:,2)=b(:,2);
x = b;
x(:,2)=0;

% x(25:end,2) = linsolve(A(25:end,25:end),b(25:end,2));
x(25:end-1,2) = -2*linsolve(A(25:end-1,25:end-1),b(25:end-1,2));

x(1:24,2)=x(25,2); % that's for the shift of beginning
x(1:89,2)=x(1:89,2)-x(90,2);
x(90:end,2)=0; % that's for the cutoff

save '%NAME.dpot.new' x '-ascii'
quit
