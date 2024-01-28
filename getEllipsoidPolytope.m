function Fk_ebox = getEllipsoidPolytope(n_iter, scaleRatio, n_facet, set_Ek, maxSegment)

%% Ellipsoid-based outer box
[X_elp, c_elp] = minve(set_Ek.V);
[~, S_elp, V_elp] = svd(X_elp);
H_elp = [V_elp'; -V_elp'];
g_elp = [sqrt(diag(S_elp)); sqrt(diag(S_elp))];
elpBox = Polyhedron('A',H_elp,'b',g_elp);
elpBox = elpBox + c_elp;
elpBox.minHRep(); elpBox.minVRep();

backup_polytope = set_Ek;

%% Volume reduction
for k = 1: n_iter
    try
        % Extract vertices
        elpBox_vtx = elpBox.V;
        
        %
        if k == 1
            % Create additional points
            n1 = size(elpBox.V,1);
            for ii = 1 : 1 : n1-1
                for jj = n1 : -1 : ii + 1
                    for mm = 2:maxSegment % 7 % 3
                        for nn = 1:mm-1
                            elpBox_vtx = [elpBox_vtx;...
                                (elpBox.V(ii,:) - elpBox.V(jj,:)).*(nn/mm) + elpBox.V(jj,:)];
                        end
                    end
                end
            end
            % Delete repeating points
            elpBox_vtx = round(elpBox_vtx,4);
            elpBox_vtx = unique(elpBox_vtx,"rows");
    
            % Exclude points inside set E
            aa = [];
            for mm = 1:size(elpBox_vtx,1)
                if ~set_Ek.contains(elpBox_vtx(mm,:)')
                    aa = [aa; elpBox_vtx(mm,:)];
                end
            end
            elpBox_vtx = aa;
        end
        
        %}

        % proceed to reduce volume
        for kk = 1 : size(elpBox_vtx,1)
            % if ~set_Ek.contains(elpBox_vtx(kk,:)')

                % create separating hyperplanes
                % temp_plane = Polyhedron('He',set_Ek.separate(elpBox_vtx(kk,:)));
                temp_point = set_Ek.project(elpBox_vtx(kk,:)').x;
                sep = temp_point - elpBox_vtx(kk,:)';
                sep(end+1) = sep'*((temp_point - elpBox_vtx(kk,:)').*scaleRatio(k,1) + elpBox_vtx(kk,:)');
                sep = sep(:)'; 
                temp_plane = Polyhedron('He',sep);
                temp_plane = temp_plane & elpBox;
                temp_plane.minHRep(); temp_plane.minVRep();

                % create separating plane-based polytope
                temp_poly1 = Polyhedron('A',-temp_plane.Ae,'b',-temp_plane.be);
                temp_poly2 = Polyhedron('A',temp_plane.Ae,'b',temp_plane.be);

                % choose the polytope that contains the set E_k
                if all(temp_poly1.contains(set_Ek.V')) %temp_poly1.contains(c_elp)
                    temp_poly = temp_poly1;
                elseif all(temp_poly2.contains(set_Ek.V'))
                    temp_poly = temp_poly2;
                else
                    error("Hyperplanes cut set Ek"); % Maybe this will never get executed, hopefully !!!
                end 

                % If the hyperplane-based polytope contains set Ek, update
                % the temporary polytope
                if all(temp_poly.contains(set_Ek.V'))
                    Fk_temp = elpBox & temp_poly;
                    Fk_temp.minHRep(); Fk_temp.minVRep();
                else
                    Fk_temp = elpBox;
                end

                % If the number of facets of the temporary polytope <
                % threshold, update ellipsoid-based outer box
                if (size(Fk_temp.getFacet(),1) < n_facet) && (all(Fk_temp.contains(set_Ek.V')))
                    elpBox = Fk_temp;
                    backup_polytope = elpBox;
                else
                    continue
                end
                
            % end
        end
    catch
        % If, for some reasons, errors occur, take the backup polytope
        elpBox = backup_polytope;
        continue
    end
end

Fk_ebox = elpBox; Fk_ebox.minHRep(); Fk_ebox.minVRep();