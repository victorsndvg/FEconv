% jacobiano de diferentes 4-node quadrilateral 
% http://mms2.ensmp.fr/ef_paris/technologie/transparents/e_Pathology.pdf

P{1} = [-1,-1];
P{2} = [ 1,-1];
P{3} = [ 1, 1];
P{4} = [-1, 1];

Pc = [1 2 3 4; 1 2 4 3; 1 3 2 4; 1 3 4 2; 1 4 2 3; 1 4 3 2];

J = @(a1, a2, a3, a4, r, s) (a4(2)-a2(2))*(a3(1)-a1(1))-(a3(2)-a1(2))*(a4(1)-a2(1))+...
                           ((a3(2)-a4(2))*(a2(1)-a1(1))-(a2(2)-a1(2))*(a3(1)-a4(1)))*r+...
                           ((a4(2)-a1(2))*(a3(1)-a2(1))-(a3(2)-a2(2))*(a4(1)-a1(1)))*s;
                   
Jac = zeros(size(Pc,1),4);
for i = 1:size(Pc,1)
  Jac(i,:) = [J(P{Pc(i,1)}, P{Pc(i,2)}, P{Pc(i,3)}, P{Pc(i,4)}, -1, -1) ...
              J(P{Pc(i,1)}, P{Pc(i,2)}, P{Pc(i,3)}, P{Pc(i,4)},  1, -1) ...
              J(P{Pc(i,1)}, P{Pc(i,2)}, P{Pc(i,3)}, P{Pc(i,4)},  1,  1) ...
              J(P{Pc(i,1)}, P{Pc(i,2)}, P{Pc(i,3)}, P{Pc(i,4)}, -1,  1)];              
end
disp([Pc Jac])