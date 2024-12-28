function handle_dry_beds(hL, hR, hmin=1e-6; α=100)
    # Smooth transition weights
    w_both_dry = exp(-α * hL) * exp(-α * hR)
    w_left_dry = exp(-α * hL) * (1 - exp(-α * hR))
    w_right_dry = (1 - exp(-α * hL)) * exp(-α * hR)
	
	println("w_both_dry = ", w_both_dry)
	println("w_left_dry = ", w_left_dry)
	println("w_right_dry = ", w_right_dry)
	
end

hmin = 0.001

hL = 0.01 - hmin
hR = 0.0 - hmin
handle_dry_beds(hL, hR, hmin)
