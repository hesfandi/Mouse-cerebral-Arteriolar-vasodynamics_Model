function [q,nodpress,H]=flow_Boas_new(H,Param,bcnod,bctyp,bcprfl,bchd,ind) 
nodsegm=4;
dishemversion = 1;
nitmax1=Param.nitmax1;
mcvcorr=Param.mcvcorr;
optw=Param.optw;
vplas=Param.vplas;
cpar=Param.cpar;
viscpar=Param.viscpar;
bifpar=Param.bifpar;
facfp=Param.facfp;
constvisc=Param.constvisc;
% phaseseparation=Param.phaseseparation;
phaseseparation=0;
consthd=0.4;
H.Edges.hd=ones(length(H.Edges.hd),1)*consthd;
bchd=ones(length(bchd),1)*consthd;
varyviscosity=Param.varyviscosity;
nnod=length(H.Nodes.X);
nseg=length(H.Edges.EndNodes);
nodpress = ones(1,nnod);
cond = ones(1, nseg);
diam=H.Edges.D;
hd=H.Edges.hd;
lseg=H.Edges.L;

nodseg = ones(nodsegm,nnod); 
nodnod = ones(nodsegm,nnod); 
nodtyp = zeros(1,nnod); 
segnodname=H.Edges.EndNodes';

for iseg=1:nseg 

        inod = segnodname(1,iseg);	%//if node names are sequential, then no search is needed
        inodd = segnodname(2,iseg);	%//if node names are sequential, then no search is needed

		nodtyp(inod)=nodtyp(inod)+1;
		nodtyp(inodd)=nodtyp(inodd)+1;

 		nodseg(nodtyp(inod),inod) = iseg;
 		nodseg(nodtyp(inodd),inodd) = iseg;


        nodnod(nodtyp(inod),inod) = inodd;
		nodnod(nodtyp(inodd),inodd) = inod;
        


end

for inod=1:nnod
 nodpress(inod) = 50.;
end

for nseg=1:nseg
 q(nseg) = 0;
end


relax = 1;
ista=H.Edges.EndNodes(:,1);
iend=H.Edges.EndNodes(:,2);
for niter=1:nitmax1
		if(mod(niter,5) == 0) 
           relax = 0.8*relax;
        end
        for iseg=1:nseg
            qold(iseg)=q(iseg);
            hdold(iseg)=hd(iseg);
                if(varyviscosity == 1)    
%                 visc = viscor(diam(iseg),hd(iseg),mcvcorr,optw,vplas,cpar,viscpar);
                visc = viscor1(diam(iseg),hd(iseg),vplas,cpar,viscpar);
                else
                visc = constvisc;%
                end
                if (ismember(iseg,ind))
                    visc=visc*0.5;
                end
%                 vscc=[vscc visc];
%                 x=facfp*(diam(iseg)^4)/lseg(iseg)/visc;
                cond(iseg) = facfp*(diam(iseg)^4)/lseg(iseg)/visc;
%                 if( niter==nitmax1 && iseg==32)
%                     cond32=cond(iseg);
%                 end
%                 if( niter==nitmax1 && iseg==27)
%                     cond27=cond(iseg);
%                 end
        end
%         if niter==nitmax1
%             f=2;
%             figure;
%             plot(cond(ArtECsInds),'*')
%             hold on
%             plot(cond(PAreECsInds),'r*')
%             hold on
%             plot(cond(PCapECsInds),'g*')
%         end
        
        nodpress=solve_Boas_qp(Param,cond,H,bcnod,bcprfl,bctyp,bchd,nodtyp,nodnod,nodseg);

        for iseg=1:nseg
			q(iseg) = (nodpress(ista(iseg)) - nodpress(iend(iseg)))*cond(iseg);
        end
        
        %//calculate segment hematocrits
        if(phaseseparation == 1)
            [nodrank,nnodfl,nodtyp,nodout,nodseg,nodnod]=putrank_Boas(H,q,nodseg,nodnod);
            if(dishemversion == 1) 
                hd=dishem_generalized_Boas(q,hd,diam,nodrank,nodout,nnodfl,nseg,bcnod,bchd,nodtyp,nodseg,bifpar);               
            else 
                %dishem();
            end
        else
            for iseg=1:nseg
                    hd(iseg) = consthd;
            end
        end
        %//compare hd and q with previous values
        maxqerr = 0.;
		maxhderr = 0.;
		errsegq = 0;
		errseghd = 0;
        qold=q; % for convergance in 1 itteration
        for iseg=1:nseg

            qchange = q(iseg) - qold(iseg);
            hdchange = hd(iseg) - hdold(iseg);
            hd(iseg) = hdold(iseg) + relax*hdchange;
            if(abs(qchange) >= maxqerr)
                maxqerr = abs(qchange);
                errsegq = iseg;
            end
            if(abs(hdchange) >= maxhderr)
                maxhderr = abs(hdchange);
                errseghd = iseg;
            end
            
        end
        
        if(maxqerr < Param.qtol && maxhderr < Param.hdtol) 
%             fprintf("Flow: %i iterations\n",niter);
            break;
        end
 end
if(phaseseparation == 1  && relax < 0.9999)
        [nodrank,nnodfl,nodtyp,nodout,nodseg,nodnod]=putrank_Boas(H,q,nodseg,nodnod);
        if(dishemversion == 1) 
            hd=dishem_generalized_Boas(q,hd,diam,nodrank,nodout,nnodfl,nseg,bcnod,bchd,nodtyp,nodseg,bifpar);    
        else
            %dishem();
        end
end
 H.Edges.hd=hd;
 