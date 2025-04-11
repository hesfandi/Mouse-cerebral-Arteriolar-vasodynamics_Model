function [nodpress]=solve_Boas_qp(Param,cond,H,bcnod,bcprfl,bctyp,bchd,nodtyp,nodnod,nodseg)

nodsegm=3;
nnod=length(H.Nodes.X);
wk=ones(nodsegm,nnod);
Adj = H.adjacency;
Adj = Adj + Adj';
Adj = full(Adj);
nseg=length(H.Edges.EndNodes(:,1));
nodpress=ones(1,nnod);
nitmax=Param.nitmax;
omega=Param.omega;
tol=Param.tol;


for inod=1:nnod
    if (nodtyp(inod) > 1)
        condsum = 0;
        for i=1:nodtyp(inod)
            iseg = abs(nodseg(i,inod));
			condsum = condsum+ cond(iseg);
			wk(i,inod) = cond(iseg);
        end
        for i=1:nodtyp(inod)
             wk(i,inod) = wk(i,inod)/condsum;
        end
    end
end

nnodbc=length(bcnod);
for inodbc=1:nnodbc
		inod = bcnod(inodbc);
        if (bctyp(inodbc)==0)
            nodpress(inod) = bcprfl(inodbc);
            nodtyp(inod) = -1;
        else 
%             x = bcprfl(inodbc) / cond(nodseg(1,inod));
			wk(1,inod) = bcprfl(inodbc) / cond(nodseg(1,inod));
        end
end


for niter=1:1*nitmax
    maxerr = 0.;
    for inod=1:nnod
        if (nodtyp(inod) == 1)
            press1 = omega*(nodpress(nodnod(1,inod)) + wk(1,inod) - nodpress(inod));
        end
        if (nodtyp(inod) >= 2)
            pcondsum = 0.;
            for i=1:nodtyp(inod)
                %x=nodpress(nodnod(i,inod));
                pcondsum = pcondsum+ wk(i,inod)*nodpress(nodnod(i,inod));
            end
            press1 = omega*(pcondsum - nodpress(inod));
        end
        if (nodtyp(inod) >= 1)
            nodpress(inod) = nodpress(inod)+ press1;
            if (abs(press1) >= maxerr)
                maxerr = abs(press1);
				errnode = inod;
                
            end
            
        end
    end
    if(maxerr < tol)
%              fprintf('solved\n')
        break;
    else 
%         fprintf('*** Warning: linear iteration not converged, maxerr = %g at node %i\n', maxerr,errnode);
    end

end