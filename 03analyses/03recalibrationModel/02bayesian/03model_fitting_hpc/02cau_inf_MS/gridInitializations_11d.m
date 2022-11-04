%This function generates a matrix of initializations from grid
%(row: number of initialization; col: number of parameters)
function inits = gridInitializations_11d(lb, ub, numS, numInits, numRepeatsEachInit)
    %----------------------------------------------------------------------
    %lb, ub            : boundaries
    %numS              : number of sections (linspace(lb,ub,numS))
    %numInits          : total number of initializations
    %numRepeatsEachInit: each initialization is repeated many times
    %----------------------------------------------------------------------
    numP          = 11;
    M_comb        = [];
    inits         = [];
    M             = NaN(numP, numS);
    %find the values in between (head and tail excluded)
    for i = 1:numP; M(i,:) = linspace(lb(i),ub(i),numS); end
    
    %find all combinations using ndgrid
    [X,Y,Z,P,Q,L,V,O,U,R,S] = ndgrid(M(1,:),M(2,:),M(3,:),M(4,:),M(5,:),M(6,:),...
        M(7,:),M(8,:),M(9,:),M(10,:),M(11,:));
    for x = 1:size(X,1)
        for y = 1:size(X,2)
            for z = 1:size(X,3)
                for p = 1:size(X,4)
                    for q = 1:size(X,5)
                        for l = 1:size(X,6)
                            for m = 1:size(X,7)
                                for n = 1:size(X,8)
                                    for o = 1:size(X,9)
                                        for r = 1:size(X,10)
                                            for s = 1:size(X, 11)
                                            M_comb =[M_comb; X(x,y,z,p,q,l,m,n,o,r,s),...
                                                Y(x,y,z,p,q,l,m,n,o,r,s),Z(x,y,z,p,q,l,m,n,o,r,s),...
                                                P(x,y,z,p,q,l,m,n,o,r,s),Q(x,y,z,p,q,l,m,n,o,r,s),...
                                                L(x,y,z,p,q,l,m,n,o,r,s),V(x,y,z,p,q,l,m,n,o,r,s),...
                                                O(x,y,z,p,q,l,m,n,o,r,s),U(x,y,z,p,q,l,m,n,o,r,s),...
                                                R(x,y,z,p,q,l,m,n,o,r,s),S(x,y,z,p,q,l,m,n,o,r,s)];
                                            end
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end 
    
    %find the total number of different combinations
    numComb         = size(M_comb,1);
    %only select part of it, can't run all of them
    selectedIndices = randperm(numComb, floor(numInits/numRepeatsEachInit));
    %repeat each one numRepeatsEachInit times
    inits           = repmat(M_comb(selectedIndices,:),[numRepeatsEachInit,1]);
    if size(inits,1) ~= numInits; disp('Error!'); end
end