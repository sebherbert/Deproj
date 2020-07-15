function hts = add_plot_curvature_k2( obj, ax )
%ADD_PLOT_CURVATURE_K2 Add the tissue plot colored with the first principal curvature.
    
    epicells = obj.epicells;
    curvs = vertcat( epicells.curvatures );
    k2_curv = curvs( :, 4 );
    
    hts = add_plot_variable( obj, k2_curv, ax );
end

