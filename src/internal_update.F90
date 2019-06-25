subroutine internal_update()
  use openarray
  use config
  use variables
  implicit none
  egb = egf;etb=et;dtb=h+etb
  et = etf;dt = h + et
  utb = utf;vtb = vtf
end subroutine internal_update

