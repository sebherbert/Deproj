function to_file( obj, file_name, include_header )
%TO_FILE Exports results to a spreadsheet file.
% We include the variable units in the header.

    if nargin < 3
        include_header = true;
    end

    T = obj.to_table;
    
    if ~include_header
       writetable( T, file_name, 'WriteVariableNames', false );
       return
    end
    
    % Add a first line with units.
    [ n_rows, n_cols ] = size( T );
    header = cell( 1, n_cols );
    
    T2 = cell2table( cell( n_rows + 1, n_cols ) );
    for c = 1 : n_cols
        units =  T.Properties.VariableUnits{ c };
        var_name = T.Properties.VariableNames{ c };
        if isempty( units )
            head = var_name;
        else
            head = sprintf( '%s_(%s)', var_name, units );
        end
        header{ 1, c } = head;
    end
    
    T2{ 1, : } = header;
    T2{ 2 : end , : } = table2cell( T( : , : ) );
    writetable( T2, file_name, 'WriteVariableNames', false );
    
end

