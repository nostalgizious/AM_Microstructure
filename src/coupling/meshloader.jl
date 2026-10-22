function load_mesh(filename::String)
    nodes = []
    elements = []

    open(filename, "r") do file
        while !eof(file)
            line = strip(readline(file))

            # Skip comment lines
            if startswith(line, "**") || isempty(line)
                continue
            end

            # Process node coordinates
            if startswith(line, "*NODE")
                # Process node coordinates
                while true
                    next_line = strip(readline(file))
                    if startswith(next_line, "*") || eof(file)
                        # element_data = parse.(Int64, split(next_line, ',')[2:end])
                        # push!(elements, element_data)
                        break
                    end
                    # Append node data
                    node_data = parse.(Float64, split(next_line, ',')[2:end])
                    push!(nodes, node_data)
                end
            end
            # Process element connectivities, no regard for which region of the mesh it is in
            if startswith(line, "*ELEMENT") != 1
                while true
                    next_line = strip(readline(file))
                    if startswith(next_line, "*") || eof(file)
                        break
                    end
                    # Append node data
                    element_data = parse.(Int64, split(next_line, ',')[1:end])
                    push!(elements, element_data)
                    
                end
            end
        end
    end
    # Reshape into arrays
    node_matrix = hcat(nodes...)
    element_data_raw = hcat(elements...)
    element_matrix = element_data_raw[2:end, :]
    elem_numbers = element_data_raw[1, :]
    return node_matrix, element_matrix, elem_numbers
end