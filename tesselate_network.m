function [XsO, XuO, connO] = tesselate_network(Xs, Xu, conn, XTr, NTr)
% Code for efficiently tesselating networks in cardinal directions.
%
% Inputs
% Xs:       d x Ns matrix of specified node coordinates
% Xu:       d x Nu matrix of unspecified node coordinates
% conn:     s x 2 matrix of edges from node i to node j
% XTr:      d x 1 vector of coordinate translations per coordinate
% NTr:      d x 1 vector of number of replications per coordinate
%
% Outputs
% XsO:      d x NsO matrix of NsO tesselated specified node positions
% XuO:      d x NuO matrix of NuO tesselated unspecified node positions
% connO:    sO x 2 matrix of connections from tesslated nodes i to j

% Parameters
d = size(Xs,1);
Ns = size(Xs,2);
Nu = size(Xu,2);
s = size(conn,1);

% Matrices of transformations
XsOProx = repmat(Xs, [1, prod(NTr)]);
XuOProx = repmat(Xu, [1, prod(NTr)]);
connOProx = repmat(conn, [prod(NTr), 1]);

% Replicate
if(d==2)
    % Extend Coordinates and Connection Indices
    for i = 1:NTr(1)
        for j = 1:NTr(2)
            XsOProx(:,[1:Ns]+((i-1)*NTr(2)+(j-1))*Ns) = Xs + XTr.*[(i-1); (j-1)];
            XuOProx(:,[1:Nu]+((i-1)*NTr(2)+(j-1))*Nu) = Xu + XTr.*[(i-1); (j-1)];
            connOProx([1:s]+((i-1)*NTr(2)+(j-1))*s,:) = conn + [Ns Nu]*(j-1 + (i-1)*NTr(2));
        end
    end
elseif(d==3)
    % Extend Coordinates and Connection Indices
    for i = 1:NTr(1)
        for j = 1:NTr(2)
            for k = 1:NTr(3)
                XsOProx(:,[1:Ns]+((i-1)*NTr(2)*NTr(3)+(j-1)*NTr(3)+(k-1))*Ns) = Xs + XTr.*[(i-1); (j-1); (k-1)];
                XuOProx(:,[1:Nu]+((i-1)*NTr(2)*NTr(3)+(j-1)*NTr(3)+(k-1))*Nu) = Xu + XTr.*[(i-1); (j-1); (k-1)];
                connOProx([1:s] +((i-1)*NTr(2)*NTr(3)+(j-1)*NTr(3)+(k-1))*s,:) = conn + [Ns Nu]*(j-1 + (i-1)*NTr(2));
            end
        end
    end
end
% Filter Redundant Connections
[C,IA,IC] = unique(XsOProx','rows','stable');
XsO = C';
XuO = XuOProx;
connO = [IC(connOProx(:,1)), connOProx(:,2)];



end