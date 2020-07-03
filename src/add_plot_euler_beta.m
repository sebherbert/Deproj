function hts = add_plot_euler_beta( ax, epicells )
%ADD_PLOT_EULER_BETA Add the tissue plot colored with the 2nd euler angle.
    
    angles = vertcat( epicells.euler_angles );
    betas = rad2deg( angles( :, 2 ) );
    
    % Wrap back to 0 - 90º.
    large_betas = betas > 90;
    betas( large_betas ) = 180 - betas( large_betas );

    hts = add_plot_variable( ax, { epicells.boundary }, betas );
end

