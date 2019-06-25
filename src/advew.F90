subroutine advew()
	use openarray
	use variables
	use config
	implicit none
	wwf=(wwb*dtb*fsm-dti2*(advw-DZF(AZB(ww)*AZB(w)))*fsm)/dtf
end subroutine
