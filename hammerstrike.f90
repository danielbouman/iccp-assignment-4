module HammerStrike

  implicit none
  private

  public rigid_hammer_strike

contains


  subroutine rigid_hammer_strike(y,hammerHeight,hammerMass,hammerVelocity,hammerForce,lastHammerVelocity,compression,deltaT, &
       g,force,hammerDone,tension,bendingPrefactor,k,b,rho,deltaX)

    real*8, intent(inout) :: y(:), hammerHeight, hammerVelocity, lastHammerVelocity, compression, force(:)
    real*8, intent(inout) :: hammerForce
    real*8, intent(in)    :: g(:), hammerMass, bendingPrefactor, rho, deltaT, deltaX, tension, k, b
    real*8                :: surface, oldforce(size(force,1)), oldcompression, comppowbp1, oldcomppowbp1, netforce, &
         oldHammerHeight, pstring, pstringmin1, profile(size(g,1))
    logical, intent(out)  :: hammerDone
    logical               :: contact(size(y,1))
    integer               :: i, j, nContactElements

    pstring = 0
    pstringmin1 = 0
    nContactElements = 0
    hammerDone = .FALSE.
    oldHammerHeight = hammerHeight
    profile = g-maxval(g)
    hammerHeight = hammerHeight + hammerVelocity*deltaT + hammerForce*deltaT**2/(2*hammerMass)
    forall (i=1:size(y,1),y(i)<hammerHeight+profile(i)-compression)
       y(i) = hammerHeight+profile(i)-compression
       contact(i) = .TRUE.
       nContactElements = nContactElements+1
    end forall
    if (all(contact.eqv..FALSE., 1).eqv..TRUE.) then
       hammerDone = .TRUE.
       return
    end if

    ! Now calculate new forces:
    oldforce = force
    forall(j=1:size(y,1),contact(j).eqv..TRUE.)
       force(j) = -(y(j-2)+4*y(j-1)+7*y(j)+4*y(j+1)+y(j+2))*(1+(y(j-1)-y(j+1))**2)
       force(j) = -force(j)*bendingPrefactor
       force(j) = force(j) - tension*(y(j-1)-2*y(j)+y(j+1))
       netforce = netforce - force(j)
       hammerVelocity = hammerVelocity + (oldforce(j)+force(j))/(hammerMass*2)
    end forall
    forall(j=1:size(y,1),contact(j).eqv..TRUE.)
       pstring = pstring + rho * deltaX * hammerVelocity
    end forall
    
    oldcompression = compression
    oldcomppowbp1 = oldcompression**b

    comppowbp1 = (netforce+hammerForce)*(hammerHeight-oldHammerHeight)/2 &
         + (-pstring**2+pstringmin1**2)/(2*rho*deltaX*nContactElements) &
         - (pstring-pstringmin1)**2/(2*hammerMass)
    comppowbp1 = comppowbp1 + comppowbp1*(b+1)/k
    compression = comppowbp1**(1/(b+1))
    
    hammerHeight = hammerHeight + oldcompression - compression
    hammerForce = netforce - k*comppowbp1/compression

  end subroutine rigid_hammer_strike
  
end module HammerStrike
