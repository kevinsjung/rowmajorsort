function rowmajorsort

    % Load input varialbles
    load unordered_list.mat
    pts = input_coordinate_final;
    
    % The first point in row major can have lowest x value OR lowest y value
    
    [~, idx1] = min(pts(1,:));
    [~, idx2] = max(pts(1,:));

    ymin = pts(2, idx1);
    ymax = pts(2, idx2);

    % Initialize variables 
    startpt = zeros(2,1);
    endpt = zeros(2,1);

    % Ordered list will be same size as unordered list
    output = zeros(2,length(pts));
    outidx = 1;
    row = 1;
    
    if ymin > ymax % rows go upwards
        % Keep moving pts out of pts to output in order
        while ~isempty(pts)
            
            % Start point is the lowest y value
            [startpt(2,1), idx] = min(pts(2,:));
            startpt(1,1) = pts(1,idx);
            
            % End point is highest x value
            [endpt(1,1), idx2] = max(pts(1,:));
            endpt(2,1) = pts(2, idx2);                      

            % Function call to update pts and output an ordered row
            [pts, inline] = ptsinline(pts, startpt, endpt,1);
            
            % Update output array
            output(:,outidx:length(inline)*row) = inline;
            outidx = length(inline)*row+1;
            row = row + 1;

        end
    elseif ymax > ymin % rows go downwards
        % Keep moving pts out of pts to output in order
        while ~isempty(pts)
            % Start point is the lowest x value
            [startpt(1,1), idx] = min(pts(1,:));
            startpt(2,1) = pts(2,idx);
            
            % End point is highest y value
            [endpt(2,1), idx2] = min(pts(2,:));
            endpt(1,1) = pts(1, idx2);
          
            % Function call to update pts and output an ordered row
            [pts, inline] = ptsinline(pts, startpt, endpt,2);

            % Update output array
            output(:,outidx:length(inline)*row) = inline;
            outidx = length(inline)*row+1;
            row = row + 1;
        end
    else  % rows are parallel to X-axis     
        % Keep moving pts out of pts to output in order
        while ~isempty(pts)
            % Start point is the lowest y value
            [startpt(2,1), idx] = min(pts(2,:));
            startpt(1,1) = pts(1,idx);

            % End point is the lowest x value
            [endpt(1,1), idx2] = max(pts(1,:));
            endpt(2,1) = pts(2, idx2);

            % Function call to update pts and output an ordered row
            [pts, inline] = ptsinline(pts, startpt, endpt,1);
            
            % Update output array
            output(:,outidx:length(inline)*row) = inline;
            outidx = length(inline)*row+1;
            row = row + 1;
        end
    end
    
    % Rename variable to match input and save
    output_coordinate_final = output;
    
    save ordered_list.mat output_coordinate_final
    
    % Bonus: Visualize the ordered list
    visualizerowmajor(output_coordinate_final)

end

function [pts, inline] = ptsinline(pts, startpt, endpt, direction)
        % Function define the line between two points and identifies other
        % points that lie on the line
        m = (startpt(2) - endpt(2)) / (startpt(1) - endpt(1));
        b = startpt(2) - m*startpt(1);
        inline = [];
        indexes = [];
        
        % Iterate through all remaining points to find other points in row
        for i = 1:length(pts)
            y = m*(pts(1,i))+b;
            % Round for rounding error
            if round(y - pts(2,i),3) == 0
                % The point lies on the line
                inline = [inline, pts(:,i)];
                indexes = [indexes, i];
            end
        end
        
        % Remove matching points out of pts
        pts(:,indexes) = [];

        % Sort matching points according to orientation of picture
        if direction == 1
            inline = sortrows(inline',2)';
        elseif direction == 2
            inline = fliplr(sortrows(inline',2)');
        end

end

function visualizerowmajor(pts)
    scatter(pts(1,:), pts(2,:));
    
    % Iterate through each element and make the index the displayed number
    for i = 1:length(pts)
        b = num2str(i); c = cellstr(b);
        dx = 0.1; dy = 0.1; % displacement so the text does not overlay the data points
        text(pts(1,i)+dx, pts(2,i)+dy, c);
        hold on;
    end
end
