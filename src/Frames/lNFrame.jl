module lNFrame

using ..SIDISXSec

# functions are not exported

function get_l′T²(var::SidisVar)
    y, Q², γ² = let v=var; v.y, v.Q², v,γ² end
    return ( 1 - y - y^2 * γ² /4 ) * Q²
end

function get_PhT²(var::SidisVar)
    Mh, y, Q², γ², zh, l_dot_Ph = let v=var; v.Mh, v.y, v.Q², v.γ², v.zh, v.l_dot_Ph end
    return - Mh^2 - y^2 * γ² / Q² * l_dot_Ph^2 + 2 * y * zh * l_dot_Ph
end

end # module
