function Fk_obox = getExtremumOuterBox(set_Ek)

% Extremum-based outer box
Fk_obox = set_Ek.outerApprox();
Fk_obox.minHRep(); Fk_obox.minVRep();

end