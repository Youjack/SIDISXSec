module lNFrame

using ..SIDISXSec

# functions are not exported

function get_l′T²(var::SidisVar)
    M, xB, y, Q² = let v=var; v.M, v.xB, v.y, v.Q² end
    γ² = get_γ²(xB, Q², M)
    return ( 1 - y - y^2 * γ² /4 ) * Q²
end

function get_PhT²(var::SidisVar)
    M, Mh, xB, y, Q², zh, cosϕh, PhT²_γN = let v=var
        v.M, v.Mh, v.xB, v.y, v.Q², v.zh, v.cosϕh, v.PhT² end
    γ²  = get_γ²(xB, Q², M)
    lT²_γN = Q² * ( 1 - y - y^2 * γ² /4 )/( y^2 * (1 + γ²) )
    q_dot_Ph = 1/γ² * ( zh * Q² - √( (1 + γ²) * Q² * ( zh^2 * Q² - γ² * (Mh^2 + PhT²_γN) ) ) )
    l_dot_Ph = - √(lT²_γN * PhT²_γN) * cosϕh + ( (1 + y * γ² /2) * q_dot_Ph + (1 - y/2) * zh * Q² )/( y * (1 + γ²) )
    return - Mh^2 - y^2 * γ² / Q² * l_dot_Ph^2 + 2 * y * zh * l_dot_Ph
end

end # module
