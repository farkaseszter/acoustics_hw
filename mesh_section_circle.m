function [mesh, nodeind, elemind] = mesh_section(mesh, r, selmode)
%MESH_SECTION Return a rectangular section of a triangular mesh
%   [SECTION, NODEIND, ELEMIND] = MESH_SECTION(MESH, LIMITS, SELMODE)
%       Returns a rectangular section of the mesh MESH. The section mesh
%       SECTION contains those elements that are located within the LIMITS
%       based on the selection mode SELMODE.
% Parameters:
%  MESH    : Triangular mesh structure
%  LIMITS  : 2x3 matrix [xmin ymin zmin; xmax ymax zmax] (can contain -Inf
%            and +Inf for no selection)
%  SELMODE : Element selection criterion. The following criteria are valid.
%            'all': An element is selected if all of its nodes are within
%                   the limints. (Default selection mode)
%            'any': An element is selected if any of its nodes are within
%                   the limits.
%            'none': An element is selected if none of its nodes are within
%                   the limits.
%            'nall': An element is selected if not all of its nodes are
%                   within the limits.
%  SECTION : Triangular mesh structure containing the desired section. The
%            returned mesh contains all the nodes of the original mesh,
%            only the element structure is changed.
%  NODEIND : The indices of the selected nodes in a column vector.
%  ELEMIND : The indices of the selected elements in a column vector.
%
% See also:
%   mesh_select


if nargin < 3
    selmode = 'all';
end

% selection of appropriate nodes and elements
expression = sprintf('(x == 0) & (y.^2 + z.^2 <= (%g+0.01)^2)', r);
[nodeind, elemind] = mesh_select(mesh, expression, selmode);
mesh.elements = mesh.elements(elemind,:);
end