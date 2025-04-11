function hd=dishem_generalized_Boas(q,hd,diam,nodrank,nodout,nnodfl,nsegfl,bcnod,bchd,nodtyp,nodseg,bifpar)
isegk = 0;

nnodbc=length(bcnod);
%boundary nodes
for bnod=1: nnodbc
		inod = bcnod(bnod);
        if (nodout(inod) >= 1)
            iseg = nodseg(1,inod);
			hd(iseg) = bchd(bnod);
			isegk=isegk+1;
        end
end

%interior nodes

for in=1: nnodfl
    inod = nodrank(in);
    nodt = nodtyp(inod);
    if(nodt >= 2)
        nout = nodout(inod);
        if (nout == nodt)
            for i=1:nodt
                hd(nodseg(i,inod)) = bchd(1); %This should not normally happen
            end
        else 
            flowsum = 0.;
            hq = 0;
            diammax = 0.;
            for i=nout+1:nodt %inflow nodes
                iseg = nodseg(i,inod);
				flowsum = flowsum +abs(q(iseg));
				hq = hq +hd(iseg)*abs(q(iseg));
                if (diam(iseg) > diammax)  %maximum of the inflows
                    diammax = diam(iseg);
                end
            end
            if(nout == 1) 
                hd(nodseg(1,inod)) = hq/flowsum;	%single outflow
            else
                for ibif=1:nout-1
                    seg1 = nodseg(ibif,inod);
					seg2 = nodseg(ibif+1,inod);
					hdd = (1. - hq/flowsum)/diammax;
					diaquot = (diam(seg1)/diam(seg2))^2;
					a = bifpar(3)*(diaquot - 1.)/(diaquot + 1.)*hdd;	
					b = 1. + bifpar(2)*hdd;
					x0 = bifpar(1)*hdd;
					flow1 = abs(q(seg1));
					qikdash = (flow1/flowsum - x0)/(1. - 2.*x0);
                    if (qikdash <= 0.)
                        rbcrat = 0.;
                    elseif (qikdash >= 1.)
                        rbcrat = 1.;
                    else 
                        rbcrat = 1./(1. + exp(-a-b*log(qikdash/(1.-qikdash))));
                    end
                    hd(seg1) = rbcrat*hq/flow1;
                    if (ibif < nout-1)
                        flowsum = flowsum- flow1;		%total flow in remaining outflows
						hq = (1. - rbcrat)*hq;	%total red cell flux in remaining outflows
						diammax = diam(seg2);	%inflow diameter to next bifurcation
                    else
                        hd(seg2) = (1.-rbcrat)*hq/abs(q(seg2));  %last remaining outflow
                    end
                end
            end
            
        end
        isegk = isegk+nout;
    end
end
if(isegk ~= nsegfl)
		fprintf("*** Error in dishem, %i of %i segments processed\n",isegk,nsegfl);
end

end