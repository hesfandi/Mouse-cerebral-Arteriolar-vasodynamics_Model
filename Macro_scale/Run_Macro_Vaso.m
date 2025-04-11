
clear 
clc
close all
%%
Tmax = 600;                % [s] Total time of the simulation
dt=0.2;
t=0:dt:Tmax-dt;
BCC=60;
%%
t_MT=5;
alpha=0.5;
K=0.003;
RC=1e-3;
%%
load('NA1.mat')
% t1=t-1.8;
% % load('th.mat')
% % t1=linspace(0,599,600);
% t2=0:0.2:600;
% gamma_int=interp1(t1, gamma, t2, 'linear');
% gamma_int=[0 gamma_int(2:end-1)];
% Tmax = 700;                % [s] Total time of the simulation
% dt=0.2;
% t=0:dt:Tmax-dt;
% MT_Inh=zeros(1,length(t));
% MT_Inh(40/dt:639.9/dt)=gamma_int;

%%
Param=input_Boas();
load('Graph_vaso_final2.mat');
PA=find(H.Edges.Type==2 | H.Edges.Type == 7);
NOC=length(PA);

Z_Edges=(H.Nodes.Z(H.Edges.EndNodes(:,1))+H.Nodes.Z(H.Edges.EndNodes(:,1)))/2;
Y_Edges=(H.Nodes.Y(H.Edges.EndNodes(:,1))+H.Nodes.Y(H.Edges.EndNodes(:,1)))/2;
X_Edges=(H.Nodes.X(H.Edges.EndNodes(:,1))+H.Nodes.X(H.Edges.EndNodes(:,1)))/2;
% [q,nodpress1,~]=flow_Boas_new(H,Param,BC(:,1),BC(:,2),[BCC(3);10],BC(:,4),[]); 
% seg=(nodpress1(:,H.Edges.EndNodes(:,1)')+nodpress1(:,H.Edges.EndNodes(:,2)'))/2;

R0=H.Edges.D(PA)*0.9;
h0=0.1*H.Edges.D(PA);

PAs=[11;162;190;218;246;274;302;330;358;386;414;442;470;498;526];

%% contractility

BC(1,3)=40;
[q,nodpress,~]=flow_Boas_new(H,Param,BC(:,1),BC(:,2),BC(:,3),BC(:,4),[]);
seg=(nodpress(:,H.Edges.EndNodes(:,1)')+nodpress(:,H.Edges.EndNodes(:,2)'))/2;
H.Edges.ctl=seg'/(seg(15));

%%
tot_cell_number = length(PA);
RR=zeros(tot_cell_number,tot_cell_number);
for i = 1:tot_cell_number
    if (i==37)
        tt=1;
    end
    ista = H.Edges.EndNodes(PA(i),1);
    iend= H.Edges.EndNodes(PA(i),2);
    edges1=find(H.Edges.EndNodes(PA(:),1)==iend);
    edges2=find(H.Edges.EndNodes(PA(:),1)==ista);
    edges3=find(H.Edges.EndNodes(PA(:),2)==ista);
    edges4=find(H.Edges.EndNodes(PA(:),2)==iend);
    
    edges=[edges1', edges2', edges3', edges4'];
    edges=unique(edges);
    for j=1:length(edges)
        if (H.Edges.Type(PA(edges(j)))==2 && H.Edges.Type(PA(i))==2)
            RR(i,edges(j))=1;
        elseif (H.Edges.Type(PA(edges(j)))==7 || H.Edges.Type(PA(i))==7)
            RR(i,edges(j))=1;
        else
        RR(i,i)=0;
        end
    end
end

neighbors = zeros(tot_cell_number, max(H.degree));
for i = 1:tot_cell_number
     num=1;
     for j = 1:tot_cell_number
         if RR(i,j) ~= 0
            neighbors(i,num)=j;
            num=num+1;
         end
     end
end

R_gj = 1 *RR;  % [Gohm]  gap junctional resistance



%%

MTT=[-0.152 ];

MT=MTT(1);
BC(1,3)=BCC(1);
H.Edges.D=H.Edges.D+MT*(H.Edges.D).*H.Edges.ctl;

[q,nodpress1,~]=flow_Boas_new(H,Param,BC(:,1),BC(:,2),[BCC(1);10],BC(:,4),[]); 
seg=(nodpress1(:,H.Edges.EndNodes(:,1)')+nodpress1(:,H.Edges.EndNodes(:,2)'))/2;

P=zeros(length(R0),length(t));
Q=zeros(length(R0),length(t));
R=zeros(length(R0),length(t));
WT=zeros(length(h0),length(t));
PD=zeros(length(R0),length(t));
VV=zeros(length(R0),length(t));
VVV=zeros(length(R0),length(t));
MT_D=zeros(length(R0),length(t));
               
R(:,1)=H.Edges.D(PA);
P(:,1)=seg(PA)';
Q(:,1)=abs(q(PA)');
WT(:,1)=R(:,1).*P(:,1)./h0;
V(:,1)=-45;

WT_th=WT(:,1);

DR=0.15;
dt2=1e-4;
%%
Q4_avg=zeros(length(BCC),1);
RK_all=[];
for j=1:length(BCC)
cnt=1;    
    for i=2:length(t)
        
            H.Edges.D(PA)=R(:,i-1);
            [q,nodpress1,~]=flow_Boas_new(H,Param,BC(:,1),BC(:,2),[BCC(j);10],BC(:,4),[]); 
            seg=(nodpress1(:,H.Edges.EndNodes(:,1)')+nodpress1(:,H.Edges.EndNodes(:,2)'))/2;
            P(:,i)=seg(PA)';
            Q(:,i)=abs(q(PA)');
            
   
            
            
            if (i>5/dt)
                %% damping
               VV(:,i)=-45+ 0.019*(WT(:,round(i-t_MT/dt))-WT_th);
               V=zeros(length(R0),length(1:200));
               V(:,1)=VV(:,i);
               for x=2:2000
                   I_gj1 = zeros(1,tot_cell_number);
                    for k = 2:tot_cell_number    
                        for z = 1:nnz(neighbors(k,:))
                            if (i>5/dt)
                                I_gj1(k) = I_gj1(k) + 1/R_gj(k,neighbors(k,z))*(V(k,x-1) - V(neighbors(k,z),x-1)); 
                            else
                                I_gj1(k) = I_gj1(k) + 1/R_gj(k,neighbors(k,z))*(WT(k,i-1) - WT(neighbors(k,z),i-1));
                            end
                        end    
                    end
            
                   V(:,x)=V(:,x-1)-(1/RC)*dt2*I_gj1';
               end
               VVV(:,i)=V(:,end);
               PD(:,i)= R0.*(1+alpha*(1-exp(-K*R(:,i-1).*P(:,i-1)./(h0))));
               MT_D(:,i)=(-45-VVV(:,i)).*H.Edges.ctl(PA);
               R_new = MT_D(:,i) + PD(:,i);
                
                % Implement the constraint on R's change rate
                delta_R_max_up = DR*dt .* R(:,i-1);
                delta_R_max_down = DR/2*dt .* R(:,i-1);

                positive_changes = (R_new - R(:,i-1)) > delta_R_max_up;
                negative_changes = (R_new - R(:,i-1)) < -delta_R_max_down;



                R_new(positive_changes) = R(positive_changes,i-1) + delta_R_max_up(positive_changes);
                R_new(negative_changes) = R(negative_changes,i-1) - delta_R_max_down(negative_changes);
                
                %%
                R(:,i) = R_new;
                
            else 
                R(:,i) = R0.*(1+alpha*(1-exp(-K*R(:,i-1).*P(:,i-1)./(h0))));
                V(:,i)=-45;        
            end
            WT(:,i)=R(:,i).*P(:,i)./h0;
    end
plot(t,R(PAs(6:15),:))

end
%%