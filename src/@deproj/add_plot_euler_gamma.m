function hts = add_plot_euler_gamma( obj, ax )
%ADD_PLOT_EULER_GAMMA Add the tissue plot colored with the 3rd euler angle.
    
    epicells = obj.epicells;
    angles = vertcat( epicells.euler_angles );
    gammas = rad2deg( angles( :, 3 ) );
    hts = add_plot_variable( obj, gammas, ax );
end

