function vector = vectorInsertAfter(vector, position, new_elements)
    %insert given data after specified position in (generalized) vector.
    %specify position 0 to place data before index 1
    
    %version 1.0 - March 6, 2018 by Walter Roberson, roberson@hushmail.com
    
    %design decision: permit multidimensional vectors not just row and column vectors
    
    %design decision: this code treats all empty vectors the same, so it is willing to
    %append tf([]) onto a row vector of char for example
    
    %design decision: do permit position 0 to indicate inserting before existing content
    
    %design decision: do not permit position -1 or 'end' or anything else to indicate appending
    
    %design decision: result is same class as vector unless vector is empty
    %in which case it is the class of the new data
    
    %design decision: it is okay if converting the class of new data might lead to saturation,
    %like uint8(300) -> 255 or uint8(-300) -> 0 or single(1E50) -> Inf
    
    %design decision: if converting new data to class of original changes shape or size of data
    %then reshape converted data to appropriate sized vector and allow that. In particular:
    %* cast() of cellstr to char will not preserve multidimensional characteristics
    %even if there is only one char per cell
    %* cast() of cellstr to char can increase numel if there is more than one char per cell
    %but if that happens then what order should the 2D characters from the 2D char array
    %be inserted into the vector? The current code uses linear order rather than row order
    
    if ~isempty(vector) && ~isvector(squeeze(vector))
        error('First argument must be empty or vector');
    end
    if ~isscalar(position) || ~isnumeric(position) || position ~= fix(position) || position < 0
        error('position must be non-negative integer');
    end
    if length(vector) < position
        error('vector is too short (%d elements) to insert AFTER position %d', length(vector), position);
    end
    if ~isempty(new_elements) && ~isvector(squeeze(new_elements))
        error('Third argument must be empty or vector');
    end
    
    %treat empty specially. isrow() and iscolumn() are both false for all empty arrays, and
    %figuring out if generalized empty arrays are compatible along the right dimensions is a headache
    if isempty(vector)
        %position must have been 0 here or would not have passed the length check
        %this case is an exception to converting the type of the new data to the existing
        %so that people can use [] for the empty vector no matter what real data type
        vector = new_elements;
    elseif isempty(new_elements)
        vector = vector;        %#ok<ASGSL> %vector does not change if you insert emptiness into it
    else
        if isrow(vector) && isrow(new_elements)
            dim = 2;
        elseif iscolumn(vector) && iscolumn(new_elements)
            dim = 1;
        elseif ~ismatrix(vector) && isscalar(new_elements)
            dim = ndims(vector);
        elseif ~ismatrix(new_elements) > 2 && isscalar(vector)
            dim = ndims(new_elements);
        elseif ~ismatrix(vector) && ~ismatrix(new_elements) && ndims(vector) == ndims(new_elements)
            dim = ndims(vector);
        else
            error('Incompatible shapes for vector and new elements')
        end
        
        try
            new_data = cast(new_elements, class(vector));  %cell to char can change number of elements and shape
        catch ME
            %design decision: hide original error message here
            error('Elements of class "%s" cannot be converted to class "%s"', class(new_elements), class(vector));
        end
        
        if ~isequal(size(new_data), size(new_elements))
            new_data = reshape(new_data, [ones(1,dim-1), numel(new_data)]);      %even if the number of elements changed
        end
        
        try
            %need a try/catch because the new elements might not be compatible datatype
            vector = cat(dim, vector(1:position), new_data, vector(position+1:end));  %this is valid even for position 0
        catch ME
            %design decision: expose original error message here
            rethrow(ME)
        end
    end