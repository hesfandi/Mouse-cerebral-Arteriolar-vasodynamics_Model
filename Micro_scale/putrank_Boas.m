function [nodrank,nnodfl,nodtyp,nodout,nodseg,nodnod]=putrank_Boas(H,q,nodseg,nodnod)


nnod=length(H.Nodes.X);
nseg=length(H.Edges.EndNodes);
ista=H.Edges.EndNodes(:,1);
iend=H.Edges.EndNodes(:,2);
nk = ones(1,nnod); 
nodtyp = zeros(1,nnod); 
nodout = zeros(1,nnod); 



nsegfl = 0;	%added TWS 2010
for iseg=1:nseg
    
    if (q(iseg) >= 0)
        nod1 = ista(iseg);
        nod2 = iend(iseg);
    else
        nod1 = iend(iseg);
        nod2 = ista(iseg);
    end
  
    nodtyp(nod1)= nodtyp(nod1)+1;
    nodseg(nodtyp(nod1),nod1) = iseg;
    nodnod(nodtyp(nod1),nod1) = nod2;
    nodout(nod1)=nodout(nod1)+1;
    nsegfl= nsegfl+1;
end
for iseg=1:nseg
    
    if (q(iseg) >= 0)
        nod1 = ista(iseg);
        nod2 = iend(iseg);
    else
        nod1 = iend(iseg);
        nod2 = ista(iseg);
    end
  
    nodtyp(nod2)= nodtyp(nod2)+1;
    nodseg(nodtyp(nod2),nod2) = iseg;
    nodnod(nodtyp(nod2),nod2) = nod1;
    
end
%assign low ranks to inflow nodes
nnodfl = 0;
for inod=1:nnod
    nk(inod) = 0;
    if(nodtyp(inod) == 1 && nodout(inod) == 1)
        nnodfl=nnodfl+1;
        nk(inod) = 1;
        nodrank(nnodfl) = inod;
    end
end

%assign increasing ranks to downstream connected nodes
flag = 1;
while(flag==1)
    flag = 0;
    for inod=1:nnod
        ff=0;
        if (nk(inod) == 0 && nodtyp(inod) > 0)
            for j=nodout(inod)+1:nodtyp(inod)
                iseg = nodseg(j,inod);
                if(inod == iend(iseg) && (nk(ista(iseg)) == 0 || q(iseg) <= 0)) 
                    ff=1;
                    break;
                end
                if(inod == ista(iseg) && (nk(iend(iseg)) == 0 || q(iseg) >= 0)) 
                    ff=1;
                    break;
                end
            end
            if (ff==0)
            nnodfl=nnodfl+1;
			nk(inod) = 1;
			nodrank(nnodfl) = inod;
			flag = 1;
            end
        end
    end
    
end

%check for unprocessed nodes--should be none
for inod=1:nnod
    if(nodtyp(inod) ~= 0 && nk(inod) == 0)
       fprintf('*** Error: unprocessed node %i in putrank\n', inod);
    end
end

end
