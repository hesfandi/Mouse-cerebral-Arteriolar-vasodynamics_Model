function [H,BC]=Graph_dynamic()
N2=4;
N1=29;
%%
offset_P=0;
step_pi=0.8;
ER=0.6;
TZ1_D=5.5;
TZ2_D=5;
TZ3_D=4.5;
L_cap=44;
std_TZ=0.02;
std_cap1=0.03;
std_cap2=0.05;
L_PA=30;
L_Pi=36;
xx=0;
% D_L=9.4+1.5;
D_L=18*(1/7)^(1/5.5);

% D_LL=9.2-14*xx;
D_U=18;
D_UU=linspace(18,18-14*xx,4);
D_LL=(D_UU)*(1/7).^(1/5.5);
% D_LL=D_UU-6;
% D_UU=18-14*xx;
D_U1=30;
D_L1=20;
aa=(1:N1-1)';
ba=(2:N1)';
X_a=repmat(0,N1,1);
Y_a=repmat(-150,N1,1);
Z_a=(linspace(0,-840,N1))';
L_a=repmat(L_PA,length(aa),1);
Type_a=repmat(2,length(aa),1);
D_a=linspace(D_U,D_L,length(aa))';
%%
av=(N1+1:2*N1-1)';
bv=(N1+2:2*N1)';
X_v=repmat(0,N1,1);
Y_v=repmat(150,N1,1);
Z_v=(linspace(0,-840,N1))';
L_v=repmat(L_PA,length(av),1);
Type_v=repmat(1,length(av),1);
D_v=linspace(D_U1,D_L1,length(av))';
%%

%%
% Creating a directed graph


g_C1 = Micro_vessel_graph();

% % Plotting the graph
% Plotting the graph
h = plot(g_C1, 'Layout', 'layered');

% Extracting the Node positions
% X_C1 = (h.XData-5)*12.5;
% Y_C1 = (h.YData-4.5)*35.85;
% X_C1=X_C1';
% Y_C1=-Y_C1';

X_C1 = (h.XData-5)*12.5;
Y_C1 = (h.YData-4.2)*38;
X_C1=X_C1';
Y_C1=-Y_C1';

% % Convert the angle to radians
% theta = 45 * (pi / 180);
% 
% % Create the rotation matrix for 45 degrees clockwise
% R = [cos(theta) sin(theta); -sin(theta) cos(theta)];
% 
% % Create a matrix of your original coordinates
% original_coordinates = [X_C1'; Y_C1'];
% 
% % Multiply the rotation matrix by the original coordinates
% rotated_coordinates = R * original_coordinates;
% 
% % Extract the rotated x and y coordinates
% X_C1 = rotated_coordinates(1, :);
% Y_C1 = rotated_coordinates(2, :);

% figure
% plot(g_C1, 'XData', X_C1, 'YData', Y_C1,'LineWidth',0.1);

%%
n_nod_C=numnodes(g_C1);
Z_C1=repmat(Z_a(4),n_nod_C,1);

edgelist = g_C1.Edges{:,1};
ac1=2*N1+edgelist(:,1);
bc1=2*N1+edgelist(:,2);
L_c=repmat(L_cap,length(edgelist(:,1)),1);
Type_c=repmat(0,length(edgelist(:,1)),1);
D_c=repmat(7,length(edgelist(:,1)),1);
X_C2=X_C1;
Y_C2=Y_C1;
Z_C2=repmat(Z_a(8),n_nod_C,1);
ac2=2*N1+edgelist(:,1)+n_nod_C;
bc2=2*N1+edgelist(:,2)+n_nod_C;

X_C3=X_C1;
Y_C3=Y_C1;
Z_C3=repmat(Z_a(12),n_nod_C,1);
ac3=2*N1+edgelist(:,1)+2*n_nod_C;
bc3=2*N1+edgelist(:,2)+2*n_nod_C;

X_C4=X_C1;
Y_C4=Y_C1;
Z_C4=repmat(Z_a(16),n_nod_C,1);
ac4=2*N1+edgelist(:,1)+3*n_nod_C;
bc4=2*N1+edgelist(:,2)+3*n_nod_C;

X_C5=X_C1;
Y_C5=Y_C1;
Z_C5=repmat(Z_a(20),n_nod_C,1);
ac5=2*N1+edgelist(:,1)+4*n_nod_C;
bc5=2*N1+edgelist(:,2)+4*n_nod_C;

X_C6=X_C1;
Y_C6=Y_C1;
Z_C6=repmat(Z_a(24),n_nod_C,1);
ac6=2*N1+edgelist(:,1)+5*n_nod_C;
bc6=2*N1+edgelist(:,2)+5*n_nod_C;

X_C7=X_C1;
Y_C7=Y_C1;
Z_C7=repmat(Z_a(29),n_nod_C,1);
ac7=2*N1+edgelist(:,1)+6*n_nod_C;
bc7=2*N1+edgelist(:,2)+6*n_nod_C;




u1=[4;8;12;16;20;24;29;3+N1;4+N1;5+N1;7+N1;8+N1;9+N1;11+N1;12+N1;13+N1;15+N1;16+N1;17+N1;19+N1;20+N1;21+N1;23+N1;24+N1;25+N1;27+N1;28+N1;29+N1];
u2=[ac1(1);ac2(1);ac3(1);ac4(1);ac5(1);ac6(1);ac7(1);...
    bc1(end-2:end);bc2(end-2:end);bc3(end-2:end);bc4(end-2:end);bc5(end-2:end);bc6(end-2:end);bc7(end-2:end)];
Type_u=zeros(length(u1),1);
L_u=20*ones(length(u1),1);
D_u=7*ones(length(u1),1);
a_t=[aa;av;ac1;ac2;ac3;ac4;ac5;ac6;ac7;u1];
b_t=[ba;bv;bc1;bc2;bc3;bc4;bc5;bc6;bc7;u2];
EndNodes=[a_t,b_t];
Type=[Type_a;Type_v;Type_c;Type_c;Type_c;Type_c;Type_c;Type_c;Type_c;Type_u];
L=[L_a;L_v;L_c;L_c;L_c;L_c;L_c;L_c;L_c;L_u];
D=[D_a;D_v;D_c;D_c;D_c;D_c;D_c;D_c;D_c;D_u];
EdgeTable=table(EndNodes,Type,L,D);
X=[X_a;X_v;X_C1;X_C2;X_C3;X_C4;X_C5;X_C6;X_C7];
Y=[Y_a;Y_v;Y_C1;Y_C2;Y_C3;Y_C4;Y_C5;Y_C6;Y_C7];
Z=[Z_a;Z_v;Z_C1;Z_C2;Z_C3;Z_C4;Z_C5;Z_C6;Z_C7];
NodeTable=table(X,Y,Z);
G=graph(EdgeTable,NodeTable);
% plotgraph(G)


Art_ind=find(G.Edges.Type==2);
bif_nod=[4,8,12];
bif_edges1=find(ismember(G.Edges.EndNodes(:,1),bif_nod'));
bif_edges2=find(ismember(G.Edges.EndNodes(:,2),bif_nod'));
bif_edges=[bif_edges1;bif_edges2];
bif_edges=unique(bif_edges);
first_Cap=bif_edges(~ismember(bif_edges,Art_ind));
G.Edges.Type(first_Cap)=10;

first_nod=[G.Edges.EndNodes(first_Cap,1);G.Edges.EndNodes(first_Cap,2)];
first_nodd=first_nod(~ismember(first_nod,bif_nod));
bif_edges1=find(ismember(G.Edges.EndNodes(:,1),first_nodd));
bif_edges2=find(ismember(G.Edges.EndNodes(:,2),first_nodd));
bif_edges=[bif_edges1;bif_edges2];
bif_edges=unique(bif_edges);
second_Cap=bif_edges(~ismember(bif_edges,first_Cap));
G.Edges.Type(second_Cap)=3;

first_nod1=[G.Edges.EndNodes(second_Cap,1);G.Edges.EndNodes(second_Cap,2)];
first_nodd=first_nod1(~ismember(first_nod1,first_nodd));
bif_edges1=find(ismember(G.Edges.EndNodes(:,1),first_nodd));
bif_edges2=find(ismember(G.Edges.EndNodes(:,2),first_nodd));
bif_edges=[bif_edges1;bif_edges2];
bif_edges=unique(bif_edges);
Third_Cap=bif_edges(~ismember(bif_edges,second_Cap));
G.Edges.Type(Third_Cap)=4;

first_nod=[G.Edges.EndNodes(Third_Cap,1);G.Edges.EndNodes(Third_Cap,2)];
first_nodd=first_nod(~ismember(first_nod,first_nodd));
bif_edges1=find(ismember(G.Edges.EndNodes(:,1),first_nodd));
bif_edges2=find(ismember(G.Edges.EndNodes(:,2),first_nodd));
bif_edges=[bif_edges1;bif_edges2];
bif_edges=unique(bif_edges);
Fourth_Cap=bif_edges(~ismember(bif_edges,Third_Cap));
G.Edges.Type(Fourth_Cap)=5;

first_nod=[G.Edges.EndNodes(Fourth_Cap,1);G.Edges.EndNodes(Fourth_Cap,2)];
first_nodd=first_nod(~ismember(first_nod,first_nodd));
bif_edges1=find(ismember(G.Edges.EndNodes(:,1),first_nodd));
bif_edges2=find(ismember(G.Edges.EndNodes(:,2),first_nodd));
bif_edges=[bif_edges1;bif_edges2];
bif_edges=unique(bif_edges);
Fifth_Cap=bif_edges(~ismember(bif_edges,Fourth_Cap));
G.Edges.Type(Fifth_Cap)=6;

first_nod=[G.Edges.EndNodes(Fifth_Cap,1);G.Edges.EndNodes(Fifth_Cap,2)];
first_nodd=first_nod(~ismember(first_nod,first_nodd));
bif_edges1=find(ismember(G.Edges.EndNodes(:,1),first_nodd));
bif_edges2=find(ismember(G.Edges.EndNodes(:,2),first_nodd));
bif_edges=[bif_edges1;bif_edges2];
bif_edges=unique(bif_edges);
sixth_Cap=bif_edges(~ismember(bif_edges,Fifth_Cap));
G.Edges.Type(sixth_Cap)=6;

% first_nod=[G.Edges.EndNodes(sixth_Cap,1);G.Edges.EndNodes(sixth_Cap,2)];
% first_nodd=first_nod(~ismember(first_nod,first_nodd));
% bif_edges1=find(ismember(G.Edges.EndNodes(:,1),first_nodd));
% bif_edges2=find(ismember(G.Edges.EndNodes(:,2),first_nodd));
% bif_edges=[bif_edges1;bif_edges2];
% bif_edges=unique(bif_edges);
% seventh_Cap=bif_edges(~ismember(bif_edges,sixth_Cap));
% G.Edges.Type(seventh_Cap)=0;
first_nod=[G.Edges.EndNodes(sixth_Cap,1);G.Edges.EndNodes(sixth_Cap,2)];
first_nodd=first_nod(~ismember(first_nod,first_nodd));
bif_edges1=find(ismember(G.Edges.EndNodes(:,1),first_nodd));
bif_edges2=find(ismember(G.Edges.EndNodes(:,2),first_nodd));
bif_edges=[bif_edges1;bif_edges2];
bif_edges=unique(bif_edges);
seventh_Cap=bif_edges(~ismember(bif_edges,sixth_Cap));
ind_t=[1,2,6,7,11,12];
dd=1:length(seventh_Cap);
G.Edges.Type(seventh_Cap(ind_t))=6;
G.Edges.Type(seventh_Cap(~ismember(dd,ind_t)))=11;

first_nod=[G.Edges.EndNodes(seventh_Cap,1);G.Edges.EndNodes(seventh_Cap,2)];
first_nodd=first_nod(~ismember(first_nod,first_nodd));
bif_edges1=find(ismember(G.Edges.EndNodes(:,1),first_nodd));
bif_edges2=find(ismember(G.Edges.EndNodes(:,2),first_nodd));
bif_edges=[bif_edges1;bif_edges2];
bif_edges=unique(bif_edges);
eigth_Cap=bif_edges(~ismember(bif_edges,seventh_Cap));
G.Edges.Type(eigth_Cap)=11;


first_nod=[G.Edges.EndNodes(eigth_Cap,1);G.Edges.EndNodes(eigth_Cap,2)];
first_nodd=first_nod(~ismember(first_nod,first_nodd));
bif_edges1=find(ismember(G.Edges.EndNodes(:,1),first_nodd));
bif_edges2=find(ismember(G.Edges.EndNodes(:,2),first_nodd));
bif_edges=[bif_edges1;bif_edges2];
bif_edges=unique(bif_edges);
nineth_Cap=bif_edges(~ismember(bif_edges,eigth_Cap));
ind_t=[1,4,7];
dd=1:length(nineth_Cap);
G.Edges.Type(nineth_Cap)=11;
G.Edges.Type(nineth_Cap(ismember(dd,ind_t)))=12;




bif_nod=[16,20,24,29];
bif_edges1=find(ismember(G.Edges.EndNodes(:,1),bif_nod'));
bif_edges2=find(ismember(G.Edges.EndNodes(:,2),bif_nod'));
bif_edges=[bif_edges1;bif_edges2];
bif_edges=unique(bif_edges);
first_Cap=bif_edges(~ismember(bif_edges,Art_ind));
G.Edges.Type(first_Cap)=10;

first_nod=[G.Edges.EndNodes(first_Cap,1);G.Edges.EndNodes(first_Cap,2)];
first_nodd=first_nod(~ismember(first_nod,bif_nod));
bif_edges1=find(ismember(G.Edges.EndNodes(:,1),first_nodd));
bif_edges2=find(ismember(G.Edges.EndNodes(:,2),first_nodd));
bif_edges=[bif_edges1;bif_edges2];
bif_edges=unique(bif_edges);
second_Cap=bif_edges(~ismember(bif_edges,first_Cap));
G.Edges.Type(second_Cap)=3;

first_nod1=[G.Edges.EndNodes(second_Cap,1);G.Edges.EndNodes(second_Cap,2)];
first_nodd=first_nod1(~ismember(first_nod1,first_nodd));
bif_edges1=find(ismember(G.Edges.EndNodes(:,1),first_nodd));
bif_edges2=find(ismember(G.Edges.EndNodes(:,2),first_nodd));
bif_edges=[bif_edges1;bif_edges2];
bif_edges=unique(bif_edges);
Third_Cap=bif_edges(~ismember(bif_edges,second_Cap));
G.Edges.Type(Third_Cap)=4;

first_nod=[G.Edges.EndNodes(Third_Cap,1);G.Edges.EndNodes(Third_Cap,2)];
first_nodd=first_nod(~ismember(first_nod,first_nodd));
bif_edges1=find(ismember(G.Edges.EndNodes(:,1),first_nodd));
bif_edges2=find(ismember(G.Edges.EndNodes(:,2),first_nodd));
bif_edges=[bif_edges1;bif_edges2];
bif_edges=unique(bif_edges);
Fourth_Cap=bif_edges(~ismember(bif_edges,Third_Cap));
G.Edges.Type(Fourth_Cap)=5;

first_nod=[G.Edges.EndNodes(Fourth_Cap,1);G.Edges.EndNodes(Fourth_Cap,2)];
first_nodd=first_nod(~ismember(first_nod,first_nodd));
bif_edges1=find(ismember(G.Edges.EndNodes(:,1),first_nodd));
bif_edges2=find(ismember(G.Edges.EndNodes(:,2),first_nodd));
bif_edges=[bif_edges1;bif_edges2];
bif_edges=unique(bif_edges);
Fifth_Cap=bif_edges(~ismember(bif_edges,Fourth_Cap));
G.Edges.Type(Fifth_Cap)=6;

first_nod=[G.Edges.EndNodes(Fifth_Cap,1);G.Edges.EndNodes(Fifth_Cap,2)];
first_nodd=first_nod(~ismember(first_nod,first_nodd));
bif_edges1=find(ismember(G.Edges.EndNodes(:,1),first_nodd));
bif_edges2=find(ismember(G.Edges.EndNodes(:,2),first_nodd));
bif_edges=[bif_edges1;bif_edges2];
bif_edges=unique(bif_edges);
sixth_Cap=bif_edges(~ismember(bif_edges,Fifth_Cap));
G.Edges.Type(sixth_Cap)=6;

% first_nod=[G.Edges.EndNodes(sixth_Cap,1);G.Edges.EndNodes(sixth_Cap,2)];
% first_nodd=first_nod(~ismember(first_nod,first_nodd));
% bif_edges1=find(ismember(G.Edges.EndNodes(:,1),first_nodd));
% bif_edges2=find(ismember(G.Edges.EndNodes(:,2),first_nodd));
% bif_edges=[bif_edges1;bif_edges2];
% bif_edges=unique(bif_edges);
% seventh_Cap=bif_edges(~ismember(bif_edges,sixth_Cap));
% G.Edges.Type(seventh_Cap)=0;
first_nod=[G.Edges.EndNodes(sixth_Cap,1);G.Edges.EndNodes(sixth_Cap,2)];
first_nodd=first_nod(~ismember(first_nod,first_nodd));
bif_edges1=find(ismember(G.Edges.EndNodes(:,1),first_nodd));
bif_edges2=find(ismember(G.Edges.EndNodes(:,2),first_nodd));
bif_edges=[bif_edges1;bif_edges2];
bif_edges=unique(bif_edges);
seventh_Cap=bif_edges(~ismember(bif_edges,sixth_Cap));
ind_t=[1,2,6,7,11,12,16,17];
dd=1:length(seventh_Cap);
G.Edges.Type(seventh_Cap(ind_t))=6;
G.Edges.Type(seventh_Cap(~ismember(dd,ind_t)))=11;

first_nod=[G.Edges.EndNodes(seventh_Cap,1);G.Edges.EndNodes(seventh_Cap,2)];
first_nodd=first_nod(~ismember(first_nod,first_nodd));
bif_edges1=find(ismember(G.Edges.EndNodes(:,1),first_nodd));
bif_edges2=find(ismember(G.Edges.EndNodes(:,2),first_nodd));
bif_edges=[bif_edges1;bif_edges2];
bif_edges=unique(bif_edges);
eigth_Cap=bif_edges(~ismember(bif_edges,seventh_Cap));
G.Edges.Type(eigth_Cap)=11;


first_nod=[G.Edges.EndNodes(eigth_Cap,1);G.Edges.EndNodes(eigth_Cap,2)];
first_nodd=first_nod(~ismember(first_nod,first_nodd));
bif_edges1=find(ismember(G.Edges.EndNodes(:,1),first_nodd));
bif_edges2=find(ismember(G.Edges.EndNodes(:,2),first_nodd));
bif_edges=[bif_edges1;bif_edges2];
bif_edges=unique(bif_edges);
nineth_Cap=bif_edges(~ismember(bif_edges,eigth_Cap));
ind_t=[1,4,7,10];
dd=1:length(nineth_Cap);
G.Edges.Type(nineth_Cap)=11;
G.Edges.Type(nineth_Cap(ismember(dd,ind_t)))=12;



TZ1=find(G.Edges.Type==3); 
TZ2=find(G.Edges.Type==4); 
TZ3=find(G.Edges.Type==5); 
Z_Edges=(G.Nodes.Z(G.Edges.EndNodes(:,1))+G.Nodes.Z(G.Edges.EndNodes(:,1)))/2;
G.Edges.D(TZ1)=normrnd(TZ1_D,std_TZ,1,length(TZ1));
G.Edges.D(TZ1)=(((Z_Edges(TZ1)+91)/(-759))*ER).*G.Edges.D(TZ1)+G.Edges.D(TZ1);
G.Edges.D(TZ2)=normrnd(TZ2_D,std_TZ,1,length(TZ2));
G.Edges.D(TZ2)=(((Z_Edges(TZ2)+91)/(-759))*ER).*G.Edges.D(TZ2)+G.Edges.D(TZ2);
G.Edges.D(TZ3)=normrnd(TZ3_D,std_TZ,1,length(TZ3));
G.Edges.D(TZ3)=(((Z_Edges(TZ3)+91)/(-759))*ER).*G.Edges.D(TZ3)+G.Edges.D(TZ3);

G.Edges.D(find(G.Edges.Type==10))=2.9;
G.Edges.D(find(G.Edges.Type==6))=normrnd(4,std_cap1,1,length(find(G.Edges.Type==6)));
G.Edges.D(find(G.Edges.Type==11))=normrnd(5.5,std_cap2,1,length(find(G.Edges.Type==11)));
G.Edges.D(find(G.Edges.Type==12))=normrnd(6.5,std_cap2,1,length(find(G.Edges.Type==12)));



G.Edges.D=G.Edges.D;
G.Edges.L(find(G.Edges.Type==10))=5;

% K=G;

% plotgraph_Conv(G)


%%

g=G;
N = numnodes(g);
s=N+1;

a1=[find(g.Edges.EndNodes(:,2)==1) find(g.Edges.EndNodes(:,1)==1)];
K = addedge(g,1,s);
b1=[find(K.Edges.EndNodes(:,2)==1) find(K.Edges.EndNodes(:,1)==1)];
a=b1(~ismember(b1,a1));
K.Edges.L(a)=0.5;
K.Edges.Type(a)=7;
K.Edges.D(a)=D_U+step_pi+offset_P;
K.Nodes.X(s)=K.Nodes.X(1)+L_Pi;
K.Nodes.Y(s)=K.Nodes.Y(1);
K.Nodes.Z(s)=K.Nodes.Z(1);
for i=1:15
    s=N+i+1;
    t=N+i;
    a1=[find(K.Edges.EndNodes(:,2)==N+i) find(K.Edges.EndNodes(:,1)==N+i)];
    K = addedge(K,s,t);
    b1=[find(K.Edges.EndNodes(:,2)==N+i) find(K.Edges.EndNodes(:,1)==N+i)];
    a=b1(~ismember(b1,a1));
    K.Edges.L(a)=L_Pi;
    K.Edges.Type(a)=7;
    K.Edges.D(a)=D_U+offset_P+step_pi*(i+1);
    K.Nodes.X(s)=K.Nodes.X(s-1)+33;
    K.Nodes.Y(s)=K.Nodes.Y(s-1); 
    K.Nodes.Z(s)=K.Nodes.Z(s-1);
end
% BC_Nodes=s;
D_intersect=D_U+offset_P+step_pi*(i+1);
art1=s;
%%


g=K;
N = numnodes(g);
s=N+1;
step=0.5;
a1=[find(g.Edges.EndNodes(:,2)==N1+1) find(g.Edges.EndNodes(:,1)==N1+1)];
K = addedge(g,N1+1,s);
b1=[find(K.Edges.EndNodes(:,2)==N1+1) find(K.Edges.EndNodes(:,1)==N1+1)];
a=b1(~ismember(b1,a1));
K.Edges.L(a)=0.5;
K.Edges.Type(a)=8;
K.Edges.D(a)=D_U1;
K.Nodes.X(s)=K.Nodes.X(N1+1)+33;
K.Nodes.Y(s)=K.Nodes.Y(N1+1);
K.Nodes.Z(s)=K.Nodes.Z(N1+1);
for i=1:13
    s=N+i+1;
    t=N+i;
    a1=[find(K.Edges.EndNodes(:,2)==N+i) find(K.Edges.EndNodes(:,1)==N+i)];
    K = addedge(K,s,t);
    b1=[find(K.Edges.EndNodes(:,2)==N+i) find(K.Edges.EndNodes(:,1)==N+i)];
    a=b1(~ismember(b1,a1));
    K.Edges.L(a)=L_Pi;
    K.Edges.Type(a)=8;
    K.Edges.D(a)=D_U1+1*i;
    K.Nodes.X(s)=K.Nodes.X(s-1)+33;
    K.Nodes.Y(s)=K.Nodes.Y(s-1); 
    K.Nodes.Z(s)=K.Nodes.Z(s-1);
end
vein1=s;

%%
% plotgraph22(K)
x_step=100;
G=K;
for j=1:3
    EndNodes=G.Edges.EndNodes;
    Type=G.Edges.Type;
    D=G.Edges.D;
    L=G.Edges.L;
    X=G.Nodes.X;
    Y=G.Nodes.Y;
    Z=G.Nodes.Z;
    X1 = numnodes(G);
    aa=(X1+1:X1+N1-1)';
    ba=(X1+2:X1+N1)';
    X_a=repmat(0+x_step*j,N1,1);
    Y_a=repmat(-150,N1,1);
    Z_a=(linspace(0,-840,N1))';
    L_a=repmat(L_PA,length(aa),1);
    Type_a=repmat(2,length(aa),1);
    D_a=linspace(D_UU(j+1),D_LL(j+1),length(aa))';
    %%
    av=(X1+N1+1:X1+2*N1-1)';
    bv=(X1+N1+2:X1+2*N1)';
    X_v=repmat(0+x_step*j,N1,1);
    Y_v=repmat(150,N1,1);
    Z_v=(linspace(0,-840,N1))';
    L_v=repmat(L_PA,length(av),1);
    Type_v=repmat(1,length(av),1);
    D_v=linspace(D_U1,D_L1,length(av))';
    %%

    %%

    

    g_C1 = Micro_vessel_graph;
%     h = plot(g_C1, 'Layout', 'layered');

% Extracting the Node positions
    X_C1 = (h.XData-5)*12.5;
    Y_C1 = (h.YData-4.2)*37;

    X_C1=x_step*j+X_C1';
    Y_C1=-Y_C1';
    n_nod_C=numnodes(g_C1);
    Z_C1=repmat(Z_a(4),n_nod_C,1);
    % figure
    % plot(g_C1, 'XData', X_C1, 'YData', Y_C1,'LineWidth',0.1);
    edgelist = g_C1.Edges{:,1};
    ac1=X1+2*N1+edgelist(:,1);
    bc1=X1+2*N1+edgelist(:,2);
    L_c=repmat(L_cap,length(edgelist(:,1)),1);
    Type_c=repmat(0,length(edgelist(:,1)),1);
    D_c=repmat(7,length(edgelist(:,1)),1);
    X_C2=X_C1;
    Y_C2=Y_C1;
    Z_C2=repmat(Z_a(8),n_nod_C,1);
    ac2=X1+2*N1+edgelist(:,1)+n_nod_C;
    bc2=X1+2*N1+edgelist(:,2)+n_nod_C;

    X_C3=X_C1;
    Y_C3=Y_C1;
    Z_C3=repmat(Z_a(12),n_nod_C,1);
    ac3=X1+2*N1+edgelist(:,1)+2*n_nod_C;
    bc3=X1+2*N1+edgelist(:,2)+2*n_nod_C;

    X_C4=X_C1;
    Y_C4=Y_C1;
    Z_C4=repmat(Z_a(16),n_nod_C,1);
    ac4=X1+2*N1+edgelist(:,1)+3*n_nod_C;
    bc4=X1+2*N1+edgelist(:,2)+3*n_nod_C;

    X_C5=X_C1;
    Y_C5=Y_C1;
    Z_C5=repmat(Z_a(20),n_nod_C,1);
    ac5=X1+2*N1+edgelist(:,1)+4*n_nod_C;
    bc5=X1+2*N1+edgelist(:,2)+4*n_nod_C;

    X_C6=X_C1;
    Y_C6=Y_C1;
    Z_C6=repmat(Z_a(24),n_nod_C,1);
    ac6=X1+2*N1+edgelist(:,1)+5*n_nod_C;
    bc6=X1+2*N1+edgelist(:,2)+5*n_nod_C;

    X_C7=X_C1;
    Y_C7=Y_C1;
    Z_C7=repmat(Z_a(29),n_nod_C,1);
    ac7=X1+2*N1+edgelist(:,1)+6*n_nod_C;
    bc7=X1+2*N1+edgelist(:,2)+6*n_nod_C;




    u1=X1+[4;8;12;16;20;24;29;3+N1;4+N1;5+N1;7+N1;8+N1;9+N1;11+N1;12+N1;13+N1;15+N1;16+N1;17+N1;19+N1;20+N1;21+N1;23+N1;24+N1;25+N1;27+N1;28+N1;29+N1];
    u2=[ac1(1);ac2(1);ac3(1);ac4(1);ac5(1);ac6(1);ac7(1);...
    bc1(end-2:end);bc2(end-2:end);bc3(end-2:end);bc4(end-2:end);bc5(end-2:end);bc6(end-2:end);bc7(end-2:end)];
    Type_u=zeros(length(u1),1);
    L_u=20*ones(length(u1),1);
    D_u=7*ones(length(u1),1);
    a_t=[aa;av;ac1;ac2;ac3;ac4;ac5;ac6;ac7;u1];
    b_t=[ba;bv;bc1;bc2;bc3;bc4;bc5;bc6;bc7;u2];
    EndNodes1=[a_t,b_t];
    Type1=[Type_a;Type_v;Type_c;Type_c;Type_c;Type_c;Type_c;Type_c;Type_c;Type_u];
    L1=[L_a;L_v;L_c;L_c;L_c;L_c;L_c;L_c;L_c;L_u];
    D1=[D_a;D_v;D_c;D_c;D_c;D_c;D_c;D_c;D_c;D_u];
    EndNodes=[EndNodes;EndNodes1];
    Type=[Type;Type1];
    L=[L;L1];
    D=[D;D1];
    EdgeTable=table(EndNodes,Type,L,D);
    X11=[X_a;X_v;X_C1;X_C2;X_C3;X_C4;X_C5;X_C6;X_C7];
    Y11=[Y_a;Y_v;Y_C1;Y_C2;Y_C3;Y_C4;Y_C5;Y_C6;Y_C7];
    Z11=[Z_a;Z_v;Z_C1;Z_C2;Z_C3;Z_C4;Z_C5;Z_C6;Z_C7];
    X=[X;X11];
    Y=[Y;Y11];
    Z=[Z;Z11];
    NodeTable=table(X,Y,Z);
    G=graph(EdgeTable,NodeTable);
    
    Art_ind=find(G.Edges.Type==2);
    bif_nod=[X1+4,X1+8,X1+12];
    bif_edges1=find(ismember(G.Edges.EndNodes(:,1),bif_nod'));
bif_edges2=find(ismember(G.Edges.EndNodes(:,2),bif_nod'));
bif_edges=[bif_edges1;bif_edges2];
bif_edges=unique(bif_edges);
first_Cap=bif_edges(~ismember(bif_edges,Art_ind));
G.Edges.Type(first_Cap)=10;

first_nod=[G.Edges.EndNodes(first_Cap,1);G.Edges.EndNodes(first_Cap,2)];
first_nodd=first_nod(~ismember(first_nod,bif_nod));
bif_edges1=find(ismember(G.Edges.EndNodes(:,1),first_nodd));
bif_edges2=find(ismember(G.Edges.EndNodes(:,2),first_nodd));
bif_edges=[bif_edges1;bif_edges2];
bif_edges=unique(bif_edges);
second_Cap=bif_edges(~ismember(bif_edges,first_Cap));
G.Edges.Type(second_Cap)=3;

first_nod1=[G.Edges.EndNodes(second_Cap,1);G.Edges.EndNodes(second_Cap,2)];
first_nodd=first_nod1(~ismember(first_nod1,first_nodd));
bif_edges1=find(ismember(G.Edges.EndNodes(:,1),first_nodd));
bif_edges2=find(ismember(G.Edges.EndNodes(:,2),first_nodd));
bif_edges=[bif_edges1;bif_edges2];
bif_edges=unique(bif_edges);
Third_Cap=bif_edges(~ismember(bif_edges,second_Cap));
G.Edges.Type(Third_Cap)=4;

first_nod=[G.Edges.EndNodes(Third_Cap,1);G.Edges.EndNodes(Third_Cap,2)];
first_nodd=first_nod(~ismember(first_nod,first_nodd));
bif_edges1=find(ismember(G.Edges.EndNodes(:,1),first_nodd));
bif_edges2=find(ismember(G.Edges.EndNodes(:,2),first_nodd));
bif_edges=[bif_edges1;bif_edges2];
bif_edges=unique(bif_edges);
Fourth_Cap=bif_edges(~ismember(bif_edges,Third_Cap));
G.Edges.Type(Fourth_Cap)=5;

first_nod=[G.Edges.EndNodes(Fourth_Cap,1);G.Edges.EndNodes(Fourth_Cap,2)];
first_nodd=first_nod(~ismember(first_nod,first_nodd));
bif_edges1=find(ismember(G.Edges.EndNodes(:,1),first_nodd));
bif_edges2=find(ismember(G.Edges.EndNodes(:,2),first_nodd));
bif_edges=[bif_edges1;bif_edges2];
bif_edges=unique(bif_edges);
Fifth_Cap=bif_edges(~ismember(bif_edges,Fourth_Cap));
G.Edges.Type(Fifth_Cap)=6;

first_nod=[G.Edges.EndNodes(Fifth_Cap,1);G.Edges.EndNodes(Fifth_Cap,2)];
first_nodd=first_nod(~ismember(first_nod,first_nodd));
bif_edges1=find(ismember(G.Edges.EndNodes(:,1),first_nodd));
bif_edges2=find(ismember(G.Edges.EndNodes(:,2),first_nodd));
bif_edges=[bif_edges1;bif_edges2];
bif_edges=unique(bif_edges);
sixth_Cap=bif_edges(~ismember(bif_edges,Fifth_Cap));
G.Edges.Type(sixth_Cap)=6;

% first_nod=[G.Edges.EndNodes(sixth_Cap,1);G.Edges.EndNodes(sixth_Cap,2)];
% first_nodd=first_nod(~ismember(first_nod,first_nodd));
% bif_edges1=find(ismember(G.Edges.EndNodes(:,1),first_nodd));
% bif_edges2=find(ismember(G.Edges.EndNodes(:,2),first_nodd));
% bif_edges=[bif_edges1;bif_edges2];
% bif_edges=unique(bif_edges);
% seventh_Cap=bif_edges(~ismember(bif_edges,sixth_Cap));
% G.Edges.Type(seventh_Cap)=0;
first_nod=[G.Edges.EndNodes(sixth_Cap,1);G.Edges.EndNodes(sixth_Cap,2)];
first_nodd=first_nod(~ismember(first_nod,first_nodd));
bif_edges1=find(ismember(G.Edges.EndNodes(:,1),first_nodd));
bif_edges2=find(ismember(G.Edges.EndNodes(:,2),first_nodd));
bif_edges=[bif_edges1;bif_edges2];
bif_edges=unique(bif_edges);
seventh_Cap=bif_edges(~ismember(bif_edges,sixth_Cap));
ind_t=[1,2,6,7,11,12];
dd=1:length(seventh_Cap);
G.Edges.Type(seventh_Cap(ind_t))=6;
G.Edges.Type(seventh_Cap(~ismember(dd,ind_t)))=11;

first_nod=[G.Edges.EndNodes(seventh_Cap,1);G.Edges.EndNodes(seventh_Cap,2)];
first_nodd=first_nod(~ismember(first_nod,first_nodd));
bif_edges1=find(ismember(G.Edges.EndNodes(:,1),first_nodd));
bif_edges2=find(ismember(G.Edges.EndNodes(:,2),first_nodd));
bif_edges=[bif_edges1;bif_edges2];
bif_edges=unique(bif_edges);
eigth_Cap=bif_edges(~ismember(bif_edges,seventh_Cap));
G.Edges.Type(eigth_Cap)=11;


first_nod=[G.Edges.EndNodes(eigth_Cap,1);G.Edges.EndNodes(eigth_Cap,2)];
first_nodd=first_nod(~ismember(first_nod,first_nodd));
bif_edges1=find(ismember(G.Edges.EndNodes(:,1),first_nodd));
bif_edges2=find(ismember(G.Edges.EndNodes(:,2),first_nodd));
bif_edges=[bif_edges1;bif_edges2];
bif_edges=unique(bif_edges);
nineth_Cap=bif_edges(~ismember(bif_edges,eigth_Cap));
ind_t=[1,4,7];
dd=1:length(nineth_Cap);
G.Edges.Type(nineth_Cap)=11;
G.Edges.Type(nineth_Cap(ismember(dd,ind_t)))=12;


bif_nod=[X1+16,X1+20,X1+24,X1+29];
   bif_edges1=find(ismember(G.Edges.EndNodes(:,1),bif_nod'));
bif_edges2=find(ismember(G.Edges.EndNodes(:,2),bif_nod'));
bif_edges=[bif_edges1;bif_edges2];
bif_edges=unique(bif_edges);
first_Cap=bif_edges(~ismember(bif_edges,Art_ind));
G.Edges.Type(first_Cap)=10;

first_nod=[G.Edges.EndNodes(first_Cap,1);G.Edges.EndNodes(first_Cap,2)];
first_nodd=first_nod(~ismember(first_nod,bif_nod));
bif_edges1=find(ismember(G.Edges.EndNodes(:,1),first_nodd));
bif_edges2=find(ismember(G.Edges.EndNodes(:,2),first_nodd));
bif_edges=[bif_edges1;bif_edges2];
bif_edges=unique(bif_edges);
second_Cap=bif_edges(~ismember(bif_edges,first_Cap));
G.Edges.Type(second_Cap)=3;

first_nod1=[G.Edges.EndNodes(second_Cap,1);G.Edges.EndNodes(second_Cap,2)];
first_nodd=first_nod1(~ismember(first_nod1,first_nodd));
bif_edges1=find(ismember(G.Edges.EndNodes(:,1),first_nodd));
bif_edges2=find(ismember(G.Edges.EndNodes(:,2),first_nodd));
bif_edges=[bif_edges1;bif_edges2];
bif_edges=unique(bif_edges);
Third_Cap=bif_edges(~ismember(bif_edges,second_Cap));
G.Edges.Type(Third_Cap)=4;

first_nod=[G.Edges.EndNodes(Third_Cap,1);G.Edges.EndNodes(Third_Cap,2)];
first_nodd=first_nod(~ismember(first_nod,first_nodd));
bif_edges1=find(ismember(G.Edges.EndNodes(:,1),first_nodd));
bif_edges2=find(ismember(G.Edges.EndNodes(:,2),first_nodd));
bif_edges=[bif_edges1;bif_edges2];
bif_edges=unique(bif_edges);
Fourth_Cap=bif_edges(~ismember(bif_edges,Third_Cap));
G.Edges.Type(Fourth_Cap)=5;

first_nod=[G.Edges.EndNodes(Fourth_Cap,1);G.Edges.EndNodes(Fourth_Cap,2)];
first_nodd=first_nod(~ismember(first_nod,first_nodd));
bif_edges1=find(ismember(G.Edges.EndNodes(:,1),first_nodd));
bif_edges2=find(ismember(G.Edges.EndNodes(:,2),first_nodd));
bif_edges=[bif_edges1;bif_edges2];
bif_edges=unique(bif_edges);
Fifth_Cap=bif_edges(~ismember(bif_edges,Fourth_Cap));
G.Edges.Type(Fifth_Cap)=6;

first_nod=[G.Edges.EndNodes(Fifth_Cap,1);G.Edges.EndNodes(Fifth_Cap,2)];
first_nodd=first_nod(~ismember(first_nod,first_nodd));
bif_edges1=find(ismember(G.Edges.EndNodes(:,1),first_nodd));
bif_edges2=find(ismember(G.Edges.EndNodes(:,2),first_nodd));
bif_edges=[bif_edges1;bif_edges2];
bif_edges=unique(bif_edges);
sixth_Cap=bif_edges(~ismember(bif_edges,Fifth_Cap));
G.Edges.Type(sixth_Cap)=6;

% first_nod=[G.Edges.EndNodes(sixth_Cap,1);G.Edges.EndNodes(sixth_Cap,2)];
% first_nodd=first_nod(~ismember(first_nod,first_nodd));
% bif_edges1=find(ismember(G.Edges.EndNodes(:,1),first_nodd));
% bif_edges2=find(ismember(G.Edges.EndNodes(:,2),first_nodd));
% bif_edges=[bif_edges1;bif_edges2];
% bif_edges=unique(bif_edges);
% seventh_Cap=bif_edges(~ismember(bif_edges,sixth_Cap));
% G.Edges.Type(seventh_Cap)=0;
first_nod=[G.Edges.EndNodes(sixth_Cap,1);G.Edges.EndNodes(sixth_Cap,2)];
first_nodd=first_nod(~ismember(first_nod,first_nodd));
bif_edges1=find(ismember(G.Edges.EndNodes(:,1),first_nodd));
bif_edges2=find(ismember(G.Edges.EndNodes(:,2),first_nodd));
bif_edges=[bif_edges1;bif_edges2];
bif_edges=unique(bif_edges);
seventh_Cap=bif_edges(~ismember(bif_edges,sixth_Cap));
ind_t=[1,2,6,7,11,12,16,17];
dd=1:length(seventh_Cap);
G.Edges.Type(seventh_Cap(ind_t))=6;
G.Edges.Type(seventh_Cap(~ismember(dd,ind_t)))=11;

first_nod=[G.Edges.EndNodes(seventh_Cap,1);G.Edges.EndNodes(seventh_Cap,2)];
first_nodd=first_nod(~ismember(first_nod,first_nodd));
bif_edges1=find(ismember(G.Edges.EndNodes(:,1),first_nodd));
bif_edges2=find(ismember(G.Edges.EndNodes(:,2),first_nodd));
bif_edges=[bif_edges1;bif_edges2];
bif_edges=unique(bif_edges);
eigth_Cap=bif_edges(~ismember(bif_edges,seventh_Cap));
G.Edges.Type(eigth_Cap)=11;


first_nod=[G.Edges.EndNodes(eigth_Cap,1);G.Edges.EndNodes(eigth_Cap,2)];
first_nodd=first_nod(~ismember(first_nod,first_nodd));
bif_edges1=find(ismember(G.Edges.EndNodes(:,1),first_nodd));
bif_edges2=find(ismember(G.Edges.EndNodes(:,2),first_nodd));
bif_edges=[bif_edges1;bif_edges2];
bif_edges=unique(bif_edges);
nineth_Cap=bif_edges(~ismember(bif_edges,eigth_Cap));
ind_t=[1,4,7,10];
dd=1:length(nineth_Cap);
G.Edges.Type(nineth_Cap)=11;
G.Edges.Type(nineth_Cap(ismember(dd,ind_t)))=12;


TZ1=find(G.Edges.Type==3); 
TZ2=find(G.Edges.Type==4); 
TZ3=find(G.Edges.Type==5); 
Z_Edges=(G.Nodes.Z(G.Edges.EndNodes(:,1))+G.Nodes.Z(G.Edges.EndNodes(:,1)))/2;
G.Edges.D(TZ1)=normrnd(TZ1_D,std_TZ,1,length(TZ1));
G.Edges.D(TZ1)=(((Z_Edges(TZ1)+91)/(-759))*ER).*G.Edges.D(TZ1)+G.Edges.D(TZ1);
G.Edges.D(TZ2)=normrnd(TZ2_D,std_TZ,1,length(TZ2));
G.Edges.D(TZ2)=(((Z_Edges(TZ2)+91)/(-759))*ER).*G.Edges.D(TZ2)+G.Edges.D(TZ2);
G.Edges.D(TZ3)=normrnd(TZ3_D,std_TZ,1,length(TZ3));
G.Edges.D(TZ3)=(((Z_Edges(TZ3)+91)/(-759))*ER).*G.Edges.D(TZ3)+G.Edges.D(TZ3);

G.Edges.D(find(G.Edges.Type==10))=2.9;
G.Edges.D(find(G.Edges.Type==10))=2.9;
G.Edges.D(find(G.Edges.Type==6))=normrnd(4,std_cap1,1,length(find(G.Edges.Type==6)));
G.Edges.D(find(G.Edges.Type==11))=normrnd(5.5,std_cap2,1,length(find(G.Edges.Type==11)));
G.Edges.D(find(G.Edges.Type==12))=normrnd(6.5,std_cap2,1,length(find(G.Edges.Type==12)));
    G.Edges.D=G.Edges.D;
    G.Edges.L(find(G.Edges.Type==10))=5;
end


%%
ind1=find(degree(G)==1);
seg_ind1=[find(ismember(G.Edges.EndNodes(:,1),ind1));find(ismember(G.Edges.EndNodes(:,2),ind1))];
seg_ind1=unique(seg_ind1);
indart=seg_ind1(G.Edges.Type(seg_ind1)==2);
indvein=seg_ind1(G.Edges.Type(seg_ind1)==1);
Nodesart=[G.Edges.EndNodes(indart,1);G.Edges.EndNodes(indart,2)];
Nodesart=unique(Nodesart);
Nodesart=Nodesart(ismember(Nodesart,ind1));
NodesVein=[G.Edges.EndNodes(indvein,1);G.Edges.EndNodes(indvein,2)];
NodesVein=unique(NodesVein);
NodesVein=NodesVein(ismember(NodesVein,ind1));

X_Nodesart=G.Nodes.X(Nodesart);
Y_Nodesart=G.Nodes.Y(Nodesart);
Y_Nodevein=G.Nodes.Y(NodesVein);
X_Nodevein=G.Nodes.X(NodesVein);

ind3=find(G.Edges.Type==7);
ind3=unique([G.Edges.EndNodes(ind3,1);G.Edges.EndNodes(ind3,2)]);
ind31=find(G.Edges.Type==8);
ind31=unique([G.Edges.EndNodes(ind31,1);G.Edges.EndNodes(ind31,2)]);
X_pial_art=G.Nodes.X(ind3);
Y_pial_art=G.Nodes.Y(ind3);
X_pial_vein=G.Nodes.X(ind31);
Y_pial_vein=G.Nodes.Y(ind31);
ind4=[];
ind5=[];
for i=1:length(X_Nodesart)
    [ d, ix ] = min( sqrt((X_pial_art-X_Nodesart(i)).^2+(Y_pial_art-Y_Nodesart(i)).^2));
    ind4=[ind4,ind3(find(G.Nodes.X(ind3)==X_pial_art(ix)& G.Nodes.Y(ind3)==Y_pial_art(ix)))];
    [ d, ix ] = min( sqrt((X_pial_vein-X_Nodevein(i)).^2+(Y_pial_vein-Y_Nodevein(i)).^2));
    ind5=[ind5,ind31(find(G.Nodes.X(ind31)==X_pial_vein(ix)& G.Nodes.Y(ind31)==Y_pial_vein(ix)))];
    
end
%%
RR=[Nodesart;NodesVein];
UU=[ind4';ind5'];
II=[RR,UU];
EndNodes=[G.Edges.EndNodes;II];
Type=[G.Edges.Type;repmat(2,length(Nodesart),1);repmat(1,length(NodesVein),1)];
D=[G.Edges.D;repmat(D_U,length(Nodesart),1);repmat(D_U1,length(NodesVein),1)];
L=[G.Edges.L;repmat(20,length(Nodesart),1);repmat(20,length(NodesVein),1)];
EdgeTable=table(EndNodes,Type,L,D);
X=G.Nodes.X;
Y=G.Nodes.Y;
Z=G.Nodes.Z;
NodeTable=table(X,Y,Z);
G=graph(EdgeTable,NodeTable);
K=G;
%%

g=K;
N = numnodes(g);
s=N+1;
step=0.5;
a1=[find(g.Edges.EndNodes(:,2)==art1);find(g.Edges.EndNodes(:,1)==art1)];
K = addedge(g,art1,s);
b1=[find(K.Edges.EndNodes(:,2)==art1);find(K.Edges.EndNodes(:,1)==art1)];
a=b1(~ismember(b1,a1));
K.Edges.L(a)=L_Pi;
K.Edges.Type(a)=7;
K.Edges.D(a)=D_intersect;
K.Nodes.X(s)=K.Nodes.X(art1)+27;
K.Nodes.Y(s)=K.Nodes.Y(art1)-13;
K.Nodes.Z(s)=K.Nodes.Z(art1);
for i=1:10
    s=N+i+1;
    t=N+i;
    a1=[find(K.Edges.EndNodes(:,2)==N+i) find(K.Edges.EndNodes(:,1)==N+i)];
    K = addedge(K,s,t);
    b1=[find(K.Edges.EndNodes(:,2)==N+i) find(K.Edges.EndNodes(:,1)==N+i)];
    a=b1(~ismember(b1,a1));
    K.Edges.L(a)=L_Pi;
    K.Edges.Type(a)=7;
    K.Edges.D(a)=D_intersect;
    K.Nodes.X(s)=K.Nodes.X(s-1)+27;
    K.Nodes.Y(s)=K.Nodes.Y(s-1)-13; 
    K.Nodes.Z(s)=K.Nodes.Z(s-1);
end
BC_Nodes=s;

g=K;
N = numnodes(g);
s=N+1;
step=0.5;
a1=[find(g.Edges.EndNodes(:,2)==vein1);find(g.Edges.EndNodes(:,1)==vein1)];
K = addedge(g,vein1,s);
b1=[find(K.Edges.EndNodes(:,2)==vein1);find(K.Edges.EndNodes(:,1)==vein1)];
a=b1(~ismember(b1,a1));
K.Edges.L(a)=L_Pi;
K.Edges.Type(a)=8;
K.Edges.D(a)=40;
K.Nodes.X(s)=K.Nodes.X(vein1)+27;
K.Nodes.Y(s)=K.Nodes.Y(vein1)-13;
K.Nodes.Z(s)=K.Nodes.Z(vein1);
for i=1:10
    s=N+i+1;
    t=N+i;
    a1=[find(K.Edges.EndNodes(:,2)==N+i) find(K.Edges.EndNodes(:,1)==N+i)];
    K = addedge(K,s,t);
    b1=[find(K.Edges.EndNodes(:,2)==N+i) find(K.Edges.EndNodes(:,1)==N+i)];
    a=b1(~ismember(b1,a1));
    K.Edges.L(a)=L_Pi;
    K.Edges.Type(a)=8;
    K.Edges.D(a)=40;
    K.Nodes.X(s)=K.Nodes.X(s-1)+27;
    K.Nodes.Y(s)=K.Nodes.Y(s-1)-13; 
    K.Nodes.Z(s)=K.Nodes.Z(s-1);
end
%%


Z_Edges=(K.Nodes.Z(K.Edges.EndNodes(:,1))+K.Nodes.Z(K.Edges.EndNodes(:,2)))/2;
K.Edges.D(find(K.Edges.Type==10))=(((Z_Edges(find(K.Edges.Type==10))+91)/(-759))*3.6)+3;


BC_Nodes=[BC_Nodes,s];
consthd=0.45;
K.Edges.hd=ones(length(K.Edges.D),1)*consthd;
K.Edges.ZEdges=Z_Edges;
H=K;
BC_Values=[45;10];
BC=[BC_Nodes',[0;0],BC_Values,[0.45;0.45]];
% H.Edges.Type(find(H.Edges.Type == 11 | H.Edges.Type == 12 | H.Edges.Type == 13 | H.Edges.Type == 14 ))=0;
% plotgraph_Conv(H)
% hold on
% P = [650,0,-450] ;   % you center point 
% L = [700,400,900] ;  % your cube dimensions 
% O = P-L/2 ;       % Get the origin of cube so that P is at center 
% plotcube(L,O,.2,[0 0 1]);   % use function plotcube 
PERCENTAGE=(sum(H.Edges.L(find(H.Edges.Type==2)))+sum(H.Edges.L(find(H.Edges.Type==3)))+sum(H.Edges.L(find(H.Edges.Type==4)))+sum(H.Edges.L(find(H.Edges.Type==5)))+sum(H.Edges.L(find(H.Edges.Type==10))))/(sum(H.Edges.L))
%%
close all;
end