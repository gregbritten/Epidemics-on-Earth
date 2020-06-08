function quad(a::Number, b::Number, c::Number)
    det_term = b^2 - 4*a*c
    if a == 0
        return print("a cannot equal zero...")
    elseif det_term < 0
        re_square = -b/(2*a)
        im_square = sqrt(-det_term) * inv(2 * a)
        re_square = round(re_square, sigdigits=3)
        im_square = round(im_square, sigdigits=3)
        return print("Complex squares are $(re_square) + $(im_square*im) and $(re_square) - $(im_square*im)")
    else
        square_1 = (-b + sqrt(det_term)) * inv(2 * a)
        square_2 = (-b - sqrt(det_term)) * inv(2 * a)
        return println("$square_1 and $square_2 are your answers!!")
    end
end
