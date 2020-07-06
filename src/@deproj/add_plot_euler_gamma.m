function hts = add_plot_euler_gamma( obj, ax )
%ADD_PLOT_EULER_GAMMA Add the tissue plot colored with the 3rd euler angle.
    
    epicells = obj.epicells;
    angles = vertcat( epicells.euler_angles );
    gammas = rad2deg( angles( :, 3 ) );
    
    neg_gammas = gammas < 0;
    gammas( neg_gammas ) = 180 + gammas( neg_gammas );
    
    hts = add_plot_variable( obj, gammas, ax );
end

