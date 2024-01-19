export renumber_poset!

"""
    renumber_poset!(poset)

Renumber the poset given by the increasing integer vector `poset`.
"""
function renumber_poset!(poset)
    #
    # Renumber the poset vector
    #
    
    # Check whether the vector is increasing

    len_poset = length(poset)
    is_increasing = true

    for i = 1:len_poset-1
        if poset[i] > poset[i+1]
            is_increasing = false
        end
    end

    if !is_increasing
        error("Poset must be increasing!")
    end

    # Renumber the poset

    current_value = poset[1]
    new_value = 1

    for i = 1:len_poset
        if poset[i] == current_value
            poset[i] = new_value
        else
            current_value = poset[i]
            new_value += 1
            poset[i] = new_value
        end
    end
end

