function hts = add_plot_curvature_k1( obj, ax )
%ADD_PLOT_CURVATURE_K1 Add the tissue plot colored with the first principal curvature.
    
    epicells = obj.epicells;
    curvs = vertcat( epicells.curvatures );
    k1_curv = curvs( :, 3 );
    
    hts = add_plot_variable( obj, k1_curv, ax );
end

