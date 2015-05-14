module HammerStrike

  implicit none
  private

  public rigid_hammer_strike

contains


  subroutine rigid_hammer_strike(y,hammerHeight,hammerMass,hammerVelocity,lastHammerVelocity,compression,deltaT,profile,force, &
       hammerDone,tension,bendingPrefactor,k,b)

    real*8, intent(inout) :: y, hammerHeight, hammerVelocity, lastHammerVelocity, compression, force(:)
    real*8, intent(in)    :: profile, hammerMass, bendingPrefactor
    real*8                :: surface, oldforce(size(force,1),size(force,2))
    logical, intent(out)  :: hammerDone
    logical               :: contact(size(y,1))
    integer               :: i

    hammerDone = .FALSE.
    hammerHeight = hammerHeight + hammerVelocity*deltaT + force*deltaT**2/(2*hammerMass)
    do i=1:size(y,1)
       y(i) = max(y(i),hammerHeight+profile(i)-compression)
       contact(i) = .TRUE.
    end do
    if all(contact.eq..FALSE., 1) then
       hammerDone = .TRUE.
       return
    end if
    

    ! Now calculate new force:
    oldforce = force
    forall(j=1:size(y,1),contact = .TRUE.)
       force(j) = (y(j-2)+4*y(j-1)+7*y(j)+4*y(j+1)+y(j+2))*(1+(y(j-1)-y(j+1))**2)
       force(j) = force(j)*bendingPrefactor
       force(j) = force(j) + tension*(y(j-1)-2y(j)+y(j+1))
       netforce = force(j)
    end forall
    momentum
    hammerHeight = hammerHeight - compression
