subroutine advctw()
	use openarray
	use variables
	use config
	implicit none
	advw=DXF(AZB(u)*AXB(ww*dt))+DYF(AZB(v)*AYB(ww*dt))
end subroutine
