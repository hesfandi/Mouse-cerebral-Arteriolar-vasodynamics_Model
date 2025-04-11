function G = Micro_vessel_graph()

G = digraph;

% Adding edges for diverging bifurcations
% Level 1
G = addedge(G, 1, 2);
% Level 2
G = addedge(G, 2, 3);
G = addedge(G, 2, 4);
% Level 2
G = addedge(G, 3, 5);
G = addedge(G, 3, 6);
G = addedge(G, 4, 7);
G = addedge(G, 4, 8);
% Level 3
G = addedge(G, 5, 9);
G = addedge(G, 5, 10);
G = addedge(G, 6, 11);
G = addedge(G, 6, 12);
G = addedge(G, 7, 13);
G = addedge(G, 7, 14);
G = addedge(G, 8, 15);
G = addedge(G, 8, 16);

% Adding edges for converging bifurcations
% Level 4
G = addedge(G, 9, 17);
G = addedge(G, 10, 18);
G = addedge(G, 11, 19);
G = addedge(G, 12, 19);
G = addedge(G, 13, 20);
G = addedge(G, 14, 20);
G = addedge(G, 15, 21);
G = addedge(G, 16, 21);
G = addedge(G, 17, 22);
G = addedge(G, 18, 22);
G = addedge(G, 19, 23);
G = addedge(G, 20, 24);
G = addedge(G, 21, 25);
G = addedge(G, 22, 26);
G = addedge(G, 23, 26);
G = addedge(G, 24, 27);
G = addedge(G, 25, 28);
end