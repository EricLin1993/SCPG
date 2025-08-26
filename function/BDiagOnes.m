function [ B ] = BDiagOnes( N1,n )
    B = zeros(N1);
    for it=0:n
         B = B + diag(ones(N1-it,1),it) + diag(ones(N1-it,1),-it)  ;
    end    

end

