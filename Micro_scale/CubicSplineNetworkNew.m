function CubicSplineNetworkNew(G,Q, KStimInd)

A = G.adjacency;
D = normalize(G.Edges.D,'range',[1 8]);
% D=D/10;

% D(G.Edges.Type == 1) = 6;
% % D(G.Edges.Type == 3) = 1;
% % D(G.Edges.Type == 2) = 2;
% D(G.Edges.Type == 3) = 4;
% % D(G.Edges.Type == 4) = 1;
% % D(G.Edges.Type == 5) = 5;
REdges = D;

npts = length(A);
nodenames = cell(1,npts);

npts2 = length(G.Edges.D);
edgenames = cell(1,npts2);

for ii = 1:npts
    nodenames{ii} = num2str(ii);
end

for ii = 1:npts2
    edgenames{ii} = num2str(ii);
end
xsl_0=25;
% xsl_1=1000;
% x_max=sqrt((xsl_1-xsl_0)^2);
% y_min=sqrt((xsl_1-xsl_0)^2);

AttributeList = ...
    {'EndNodes' 'Weight' 'REdges', 'Q'};
edgenod = G.Edges{:,1};
nseg = numedges(G);
EdgeTable = table([edgenod(:,1) edgenod(:,2)],ones(nseg,1),REdges, Q,...
    'VariableNames',AttributeList);

X = G.Nodes.X; Y = G.Nodes.Y; Z = G.Nodes.Z;
NodeTable = table(X,Y,Z,nodenames', 'VariableNames',{'X' 'Y' 'Z' 'nodenames'});
g = graph(EdgeTable,NodeTable);
xyz = [g.Nodes.X, g.Nodes.Y, g.Nodes.Z];


%% remove edges of the nodes with degree higher than 2
% find all nodes with degrees higher than 2
Bif_nodeIDs = find(g.degree>2);
g_modified = g;
for ii = 1:numel(Bif_nodeIDs)
    g_modified = g_modified.rmedge(g_modified.outedges(Bif_nodeIDs(ii)));
end

%%
[bins,binsizes] = g_modified.conncomp('outputform','cell');
strands = bins;

%%

remove_indx = [];
badstrands = [];

for ii = 1:numel(strands)
    strand = strands{ii};
    if isequal(numel(strand), 1)
        continue
    end
    
    End_nodes = strand((ismember(g_modified.degree(strand),[0,1])));
    full_strand = strand;
    if ~isempty(End_nodes)
        for jj = 1:numel(End_nodes)
            Add_Neighbor = intersect(g.neighbors(End_nodes(jj)), Bif_nodeIDs);
            full_strand = [Add_Neighbor', full_strand];
        end
    end
    full_strand = unique(full_strand,'stable');
    strands{ii} = full_strand;
    
    % check if the order is okay
    check = 1;
    istrand = strands{ii};
    for i = 1:numel(istrand) - 1
        nod1 = istrand(i);
        nod2 = istrand(i + 1);
        iedge = findedge(g, nod1, nod2);
        if isequal(iedge, 0)
            check = 0;
            break
        end
    end
    
    if check
        continue
    end
    
    sg = subgraph(g,full_strand);
    
    % find degree ones of the strand
    NodeDegrees = sg.degree(1:numel(full_strand));
    Deg1Nodes = find(NodeDegrees == 1);
    if isempty(Deg1Nodes)
        % loop is happening
        % remove one egde
        nodeIDs = sg.Edges{1,1};
        node1 = nodeIDs(1);
        node2 = nodeIDs(2);
        sg1 = sg;
        sg1 = rmedge(sg1,node1,node2);
        full_strand = full_strand(shortestpath(sg1,node1,node2));
        strands{ii} = full_strand;
        continue
    end
    startNode = Deg1Nodes(1);
    new_strand = [startNode];
    sg_modified = sg;
    
    while length(new_strand)<length(full_strand)
        % find neighbors of start node
        neighbor_nodes = sg_modified.neighbors(startNode);
        new_strand = [new_strand, neighbor_nodes];
        sg_modified = sg_modified.rmedge(sg_modified.outedges(startNode));
        startNode = new_strand(end);
    end
    new_strand = unique(new_strand,'stable');
    strands{ii} = full_strand(new_strand);
    
end

%%
ind = find(binsizes == 1);
count = numel(strands);
for i = 1:numel(ind)
    istrand = strands{ind(i)};
    if isstring(istrand)
        strand = str2double(istrand);
    else
        strand = istrand;
    end
    Add_Neighbor = g.neighbors(strand);
    for j = 1:numel(Add_Neighbor)
        count = count + 1;
        strands{count} = [strand, Add_Neighbor(j)];
    end
end

%% plot with cubic spline

disp('plotting.......')
tic
CubicSplineFig = figure;
CubicSplineFig.WindowState = 'maximized';
% set(gca, 'ZDir','reverse')
set(gcf, 'visible','off');
hold all

n = 30;  % number of faces for tubes
inner_points = 3;    % number of points in the axial direction for each segment

errcount = 0;
errs = [];
count = 1;

for s = 1:length(strands)

    strand = strands{s};
    if numel(strand) == 1
        continue
    end
    
    XYZ = xyz(strand,:)';
    cs = cscvn(XYZ);
    for ii = 1:length(cs.breaks) - 1
        nod1 = strand(ii);
        nod2 = strand(ii + 1);
        iedge = findedge(g, nod1, nod2);
        if length(iedge) == 2
            iedge=iedge(1);
        end
        
        if iedge == 0
            errcount = errcount + 1;
            errs = [errs, s];
            continue
        else
            rval = g.Edges.REdges(iedge);
            strandQ = g.Edges.Q(iedge);
        end
        xtest = linspace(cs.breaks(ii),cs.breaks(ii+1),inner_points);
        ytest = linspace(cs.breaks(ii),cs.breaks(ii+1),inner_points);
        ztest = linspace(cs.breaks(ii),cs.breaks(ii+1),inner_points);
        
        ctestx = 0;
        ctesty = 0;
        ctestz = 0;
        
        for jj = 1:cs.order
            ctestx = ctestx + cs.coefs((ii-1)*3+1,jj)*(xtest - cs.breaks(ii)).^(cs.order - jj);
            ctesty = ctesty + cs.coefs((ii-1)*3+2,jj)*(ytest - cs.breaks(ii)).^(cs.order - jj);
            ctestz = ctestz + cs.coefs((ii-1)*3+3,jj)*(ztest - cs.breaks(ii)).^(cs.order - jj);
        end
        
        [X,Y,Z] = tubeplot([ctestx;ctesty;ctestz],rval,n,0);
        h = surface(X,Y,Z,'edgecolor','none');
        all_hs(count) = h;
        allEdgeIDs(count) = iedge;
        count = count + 1;
        
        
    end
end

fig = figure(CubicSplineFig);
fig.Color = 'w';

% shading interp
view(44,6)
brighten(1)


map = jet;
colormap(map)
% caxis([0 5])
caxis([min(g.Edges.Q) max(g.Edges.Q)])

ax = gca;
ax.Color = 'k';
ax.LineWidth = 1;
ax.FontSize = 16;
ax.FontName = 'arial';
ax.FontWeight = 'bold';
ax.XTick = [];
ax.YTick = [];
ax.ZTick = [];
axis image;
axis equal
toc

%% assign color for edges
for i = 1:numel(all_hs)
    h = all_hs(i);
    iedge = allEdgeIDs(i);
    strandQ = g.Edges.Q(iedge);
    h.CData = strandQ*ones(size(h.ZData));
end

%% add spheres at stimulated nodes
stimNodes = KStimInd;
for i = 1:length(stimNodes)
    nodeID = stimNodes(i);
    rs = 3*mean(g.Edges.REdges);          % radius for the sphere representing the cells
    
    [Xs,Ys,Zs] = sphere(10);
    
    Xs = rs*Xs + g.Nodes.X(nodeID);
    Ys = rs*Ys + g.Nodes.Y(nodeID);
    Zs = rs*Zs + g.Nodes.Z(nodeID);
    
    h = surf(Xs,Ys,Zs,'linestyle','none');
    h.FaceColor = 'w';
end

%% add text( Edge number) at Edges
% % for i = 1:length(G.Edges.D)
% %     nodeID1 = G.Edges.EndNodes(i,1);
% %     nodeID2 = G.Edges.EndNodes(i,2);
% %     
% %     Xs = ((G.Nodes.X(nodeID1)+G.Nodes.X(nodeID2))/2);
% %     Ys = ((G.Nodes.Y(nodeID1)+G.Nodes.Y(nodeID2))/2);
% %     Zs = (G.Nodes.Z(nodeID1)+G.Nodes.Z(nodeID2))/2;
% %     
% %     text(Xs,Ys,Zs,edgenames{i},'Color','w','FontSize',3,'FontWeight','bold');
% % end
% view(2)
colorbar
end

function [x,y,z]=tubeplot(curve,r,n,ct)
% Usage: [x,y,z]=tubeplot(curve,r,n,ct)
%
% Tubeplot constructs a tube, or warped cylinder, along
% any 3D curve, much like the build in cylinder function.
% If no output are requested, the tube is plotted.
% Otherwise, you can plot by using surf(x,y,z);
%
% Example of use:
% t=linspace(0,2*pi,50);
% tubeplot([cos(t);sin(t);0.2*(t-pi).^2],0.1);
% daspect([1,1,1]); camlight;
%
% Arguments:
% curve: [3,N] vector of curve data
% r      the radius of the tube
% n      number of points to use on circumference. Defaults to 8
% ct     threshold for collapsing points. Defaults to r/2
%
% The algorithms fails if you have bends beyond 90 degrees.
% Janus H. Wesenberg, july 2004

if nargin<3 || isempty(n), n=8;
    if nargin<2, error('Give at least curve and radius');
    end;
end;
if size(curve,1)~=3
    error('Malformed curve: should be [3,N]');
end;
if nargin<4 || isempty(ct)
    ct=0.5*r;
end


%Collapse points within 0.5 r of each other
npoints=1;
for k=2:(size(curve,2)-1)
    if norm(curve(:,k)-curve(:,npoints))>ct;
        npoints=npoints+1;
        curve(:,npoints)=curve(:,k);
    end
end
%Always include endpoint
if norm(curve(:,end)-curve(:,npoints))>0
    npoints=npoints+1;
    curve(:,npoints)=curve(:,end);
end

%deltavecs: average for internal points.
%           first strecth for endpoitns.
dv=curve(:,[2:end,end])-curve(:,[1,1:end-1]);

%make nvec not parallel to dv(:,1)
nvec=zeros(3,1);
[buf,idx]=min(abs(dv(:,1))); nvec(idx)=1;

xyz=repmat([0],[3,n+1,npoints+2]);

%precalculate cos and sing factors:
cfact=repmat(cos(linspace(0,2*pi,n+1)),[3,1]);
sfact=repmat(sin(linspace(0,2*pi,n+1)),[3,1]);

%Main loop: propagate the normal (nvec) along the tube
for k=1:npoints
    convec=cross(nvec,dv(:,k));
    convec=convec./norm(convec);
    nvec=cross(dv(:,k),convec);
    nvec=nvec./norm(nvec);
    %update xyz:
    xyz(:,:,k+1)=repmat(curve(:,k),[1,n+1])+...
        cfact.*repmat(r*nvec,[1,n+1])...
        +sfact.*repmat(r*convec,[1,n+1]);
end;

%finally, cap the ends:
xyz(:,:,1)=repmat(curve(:,1),[1,n+1]);
xyz(:,:,end)=repmat(curve(:,end),[1,n+1]);

%,extract results:
x=squeeze(xyz(1,:,:));
y=squeeze(xyz(2,:,:));
z=squeeze(xyz(3,:,:));

%... and plot:
if nargout<3, surf(x,y,z); end;
end