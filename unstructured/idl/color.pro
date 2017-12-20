function color, c, maxcolors
    col = get_colors(maxcolors)
    return, col[c mod n_elements(col)] 
end
