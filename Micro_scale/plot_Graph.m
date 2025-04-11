function plot_Graph(q,H)
KStimInd = [];
D=H.Edges.diam;
L=H.Edges.len;
Types=H.Edges.types;
for i=1:length(Types)
    if (D(i)>=12)
        type(i)=2;
    elseif (D(i)>3 && D(i)<12)
        type(i)=1;
    else
        type(i)=0;
    end
end
Type=type';
X=H.Nodes.x;
Y=H.Nodes.y;
Z=H.Nodes.z;
EndNodes=H.Edges.EndNodes;
EdgeTable=table(EndNodes,D,L,Type);
NodeTable=table(X,Y,Z);
T=graph(EdgeTable,NodeTable);
climits = [0 500];
Q=q';
CubicSplineNetworkNew(T,Q, KStimInd)
set(gcf, 'InvertHardcopy', 'off')
end