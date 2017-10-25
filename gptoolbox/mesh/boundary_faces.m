function [F,J,K] = boundary_faces(T)
  % BOUNDARY_FACES Determine boundary faces of tetrahedra stored in T
  %
  % F = boundary_faces(T)
  %
  % Input:
  %   T  tetrahedron index list, m by 4, where m is the number of tetrahedra
  % Output:
  %   F  list of boundary faces, n by 3, where n is the number of boundary faces
  %   J  list of indices into T, n by 1
  %   K  list of indices revealing across from which vertex is this face
  %

  % get all faces
  allF = [ ...
    T(:,[4 2 3]); ...
    T(:,[3 1 4]); ...
    T(:,[2 4 1]); ...
    T(:,[1 3 2])];
  % sort rows so that faces are reorder in ascending order of indices
  sortedF = sort(allF,2);

  % This is a wild hack but shaves nearly a 2x speed up. Convert the rows to
  % integers before calling unique
  if all(sortedF(:))>=0
    sortedF = uint64(sortedF);
    max_f = max(sortedF(:))+1;
    max_value = intmax('uint64');
    if max_f*max_f*max_f < max_value;
      sortedF = [sortedF(:,1)+sortedF(:,2)*max_f+sortedF(:,3)*max_f^2];
    elseif max_f*max_f < max_value
      sortedF = [sortedF(:,1)+sortedF(:,2)*max_f sortedF(:,3)];
    end
  end

  % determine uniqueness of faces
  [u,m,n] = unique(sortedF,'rows');
  % determine counts for each unique face
  counts = accumarray(n(:), 1);
  % extract faces that only occurred once
  I = m(counts == 1);
  L = I-1;
  J = mod(L,size(T,1))+1;
  K = floor(L/size(T,1))+1;
  F = allF(I,:);
end
