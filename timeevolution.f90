module TimeEvolution

implicit none
public rigid_hammer_strike

public time_evolution
public time_evolution_bridge

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine time_evolution(y,initY,g,N,nTimeSteps,Ms, &
 bL1,bL2,bL3,bL4,bLF,a1,a2,a3,a4,a5,bR1,bR2,bR3,bR4,bRF,deltaT)


  real*8, intent(in)  :: initY(:), g(:), Ms, deltaT

  real*8, intent(in)  :: bL1,bL2,bL3,bL4,bLF
  real*8, intent(in)  :: a1,a2,a3,a4,a5
  real*8, intent(in)  :: bR1,bR2,bR3,bR4,bRF
  integer, intent(in) :: N, nTimeSteps

  real*8, intent(out) :: y(N+1,nTimeSteps+1)

  real*8              :: F
  integer             :: i, t


  !!! First step 1: !!!
  t = 1
  F = 0

    y(1,1) = bL1*initY(1)+bL2*initY(2)+bL3*initY(3)+bL4*initY(1)+bLF*0
    y(2,1) = a1*(initY(4)-initY(2)+2*initY(1)) + a2*(initY(3)+initY(1)) &
 + a3*initY(2) + a4*initY(2) + a5*(initY(3)+initY(1))
    y(N-1,1) = a1*(2*initY(N)-initY(N-1)+initY(N-3)) &
 + a2*(initY(N)+initY(N-2)) + a3*initY(N-1) + a4*initY(N-1) &
 + a5*(initY(N)+initY(N-2))
    y(N,1) = bR1*initY(N) + bR2*initY(N-1) + bR3*initY(N-1) + bR4*initY(N) &
 + bRF*0

    do i=3,(N-1)
      y(i,t) = a1*(initY(i+2)+initY(i-2)) + a2*(initY(i+1) + initY(i-1)) &
 + a3*initY(i) + a4*initY(i) + a5*(initY(i+1)+initY(i-1)) &
 + deltaT*deltaT*N*F*g(i)/Ms
    end do

 ! Now the second step:
 t = 2
   
    y(1,t) = bL1*y(1,t-1)+bL2*y(2,t-1)+bL3*y(3,t-1)+bL4*initY(1)+bLF*0
    y(2,t) = a1*(y(4,t-1)-y(2,t-1)+2*y(1,t-1)) + a2*(y(3,t-1)+y(1,t-1)) &
 + a3*y(2,t-1) + a4*initY(2) + a5*(initY(3)+initY(1))
    y(N-1,t) = a1*(2*y(N,t-1)-y(N-1,t-1)+y(N-3,t-1)) &
 + a2*(y(N,t-1)+y(N-2,t-1)) + a3*y(N-1,t-1) + a4*initY(N-1) &
 + a5*(initY(N)+initY(N-2))
    y(N,t) = bR1*y(N,t-1) + br2*y(N-1,t-1) + bR3*y(N-1,t-1) + bR4*initY(N) &
 + bRF*0

    do i=3,(N-1)
      y(i,t) = a1*(y(i+2,t-1)+y(i-2,t-1)) + a2*(y(i+1,t-1) + y(i-1,t-1)) &
 + a3*y(i,t-1) + a4*initY(i) + a5*(initY(i+1)+initY(i-1)) &
 + deltaT*deltaT*N*F*g(i)/Ms
    end do


! Finally the Great Time Loop:

do t=3,(nTimeSteps+1)
    y(1,t) = bL1*y(1,t-1)+bL2*y(2,t-1)+bL3*y(3,t-1)+bL4*y(1,t-2)+bLF*0
    y(2,t) = a1*(y(4,t-1)-y(2,t-1)+2*y(1,t-1)) + a2*(y(3,t-1)+y(1,t-1)) &
 + a3*y(2,t-1) + a4*y(2,t-2) + a5*(y(3,t-2)+y(1,t-2))
    y(N-1,t) = a1*(2*y(N,t-1)-y(N-1,t-1)+y(N-3,t-1)) &
 + a2*(y(N,t-1)+y(N-2,t-1)) + a3*y(N-1,t-1) + a4*y(N-1,t-2) &
 + a5*(y(N,t-2)+y(N-2,t-2))
    y(N,t) = bR1*y(N,t-1) + br2*y(N-1,t-1) + bR3*y(N-1,t-1) + bR4*y(N,t-2) &
 + bRF*0

    do i=3,(N-1)
      y(i,t) = a1*(y(i+2,t-1)+y(i-2,t-1)) + a2*(y(i+1,t-1) + y(i-1,t-1)) &
 + a3*y(i,t-1) + a4*y(i,t-2) + a5*(y(i+1,t-2)+y(i-1,t-2)) &
 + deltaT*deltaT*N*F*g(i)/Ms
    end do

  if (t==1) then
    F = 1.5d1
  else if (t==100) then
    F = 0d0
  end if

end do

end subroutine


!!!!!!!!!!!!

subroutine time_evolution_bridge(bridgeY,initY,g,N,nTimeSteps,Ms,bridgePosition, &
     deltaT,inputF,fduration,b1,b2,deltaX,r,mu,zetaB,zetaL,rho,b1Damped,b2Damped,dampT, &
     initialHammerheight,initialHammerVelocity,k,b,hammerMass)

  real*8, intent(in)  :: initY(:), g(:), Ms, deltaT, bridgePosition, inputF, fduration, k, b
  real*8, intent(in)  :: b1, b2, deltaX, r, mu, zetaB, zetaL, rho, b1Damped, b2Damped, dampT, initialHammerHeight, &
       initialHammerVelocity, hammerMass

  real*8              :: bL1,bL2,bL3,bL4,bLF
  real*8              :: a1,a2,a3,a4,a5
  real*8              :: bR1,bR2,bR3,bR4,bRF

  real*8              :: hammerHeight(size(initY,1)), hammerVelocity(size(initY,1)), hammerForce(size(initY,1)), &
       lastHammerVelocity, compression, &
                              stringForce(size(initY,1)), bendingPrefactor, tension
  integer, intent(in) :: N, nTimeSteps
  logical             :: hammerDone

  real*8, intent(out) :: bridgeY(3,nTimeSteps+1)

  real*8              :: F, y(N+1,3), vy(N+1)
  integer             :: i, t, localt, bridgeElement

  real*8              :: D, nu, b_R_denom, b_L_denom

  y = 0
  hammerHeight = initialHammerHeight
  hammerVelocity = initialHammerVelocity
  hammerForce = 0
  stringForce = 0
  hammerDone = .FALSE.

  !!! Constants

  D = 1+b1*deltaT
  nu = (2*b2*deltaT)/(deltaX**2)
    
  a1 = (-(r**2)*mu)/D
  a2 = (r**2+4*(r**2)*mu+nu)/D
  a3 = (2-2*(r**2)-6*(r**2)*mu-2*nu)/D
  a4 = (-1+b1*deltaT+2*nu)/D
  a5 = -nu/D

  b_R_denom = 1 + b1*deltaT+zetaB*r
  bR1 = (2-2*r**2*mu - 2*r**2)/b_R_denom
  bR2 = (4*r**2*mu - 2*r**2)/b_R_denom
  bR3 = (-2*r**2*mu)/b_R_denom
  bR4 = (-1+b1*deltaT + zetaB*r)/b_R_denom
  bRF = (deltaT**2/rho)/b_R_denom
  
  b_L_denom = 1 + b1*deltaT+zetaL*r
  bL1 = (2-2*r**2*mu - 2*r**2)/b_L_denom
  bL2 = (4*r**2*mu - 2*r**2)/b_L_denom
  bL3 = (-2*r**2*mu)/b_L_denom
  bL4 = (-1+b1*deltaT + zetaL*r)/b_L_denom
  bLF = (deltaT**2/rho)/b_L_denom

  tension = (k*deltaX*rho/deltaT)**2
  bendingPrefactor = mu*tension*deltaX*deltaX
  
  bridgeElement = nint(bridgePosition*N)
  !!! First step 1: !!!
  t = 1
  localt = 3
!!!!!!!!!!!!!!!!! Magic numbers are highly discouraged, hence
!!!!!!!!!!!!!!!!! localt having a constant value of 3 (which
!!!!!!!!!!!!!!!!! helps continuity as well.)
  F = 0
  vy = 0

    y(1,localt) = bL1*initY(1)+bL2*initY(2)+bL3*initY(3)+bL4*initY(1)+bLF*0
    y(2,localt) = a1*(initY(4)-initY(2)+2*initY(1)) + a2*(initY(3)+initY(1)) &
 + a3*initY(2) + a4*initY(2) + a5*(initY(3)+initY(1))
    y(N-1,localt) = a1*(2*initY(N)-initY(N-1)+initY(N-3)) &
 + a2*(initY(N)+initY(N-2)) + a3*initY(N-1) + a4*initY(N-1) &
 + a5*(initY(N)+initY(N-2))
    y(N,localt) = bR1*initY(N) + bR2*initY(N-1) + bR3*initY(N-1) + bR4*initY(N) &
 + bRF*0

    do i=3,(N-1)
      y(i,localt) = a1*(initY(i+2)+initY(i-2)) + a2*(initY(i+1) + initY(i-1)) &
 + a3*initY(i) + a4*initY(i) + a5*(initY(i+1)+initY(i-1)) &
 + deltaT*deltaT*N*stringForce(i)/Ms
    end do

! Transfer the bridge element string amplitudes to bridgeY:
  do i=1,3
    bridgeY(i,t) = y(bridgeElement-2+i,localt)
  end do

  y(:,localt-2) = y(:,localt-1)
  y(:,localt-1) = y(:,localt)

 ! Now the second step:
 t = 2
   
    y(1,localt) = bL1*y(1,localt-1)+bL2*y(2,localt-1)+bL3*y(3,localt-1)+bL4*initY(1)+bLF*0
    y(2,localt) = a1*(y(4,localt-1)-y(2,localt-1)+2*y(1,localt-1)) + a2*(y(3,localt-1)+y(1,localt-1)) &
 + a3*y(2,localt-1) + a4*initY(2) + a5*(initY(3)+initY(1))
    y(N-1,localt) = a1*(2*y(N,localt-1)-y(N-1,localt-1)+y(N-3,localt-1)) &
 + a2*(y(N,localt-1)+y(N-2,localt-1)) + a3*y(N-1,localt-1) + a4*initY(N-1) &
 + a5*(initY(N)+initY(N-2))
    y(N,localt) = bR1*y(N,localt-1) + br2*y(N-1,localt-1) + bR3*y(N-1,localt-1) + bR4*initY(N) &
 + bRF*0

    do i=3,(N-1)
      y(i,localt) = a1*(y(i+2,localt-1)+y(i-2,localt-1)) + a2*(y(i+1,localt-1) + y(i-1,localt-1)) &
 + a3*y(i,localt-1) + a4*initY(i) + a5*(initY(i+1)+initY(i-1)) &
 + deltaT*deltaT*N*stringForce(i)/Ms
    end do

! Transfer the bridge element string amplitudes to bridgeY:
  do i=1,3
    bridgeY(i,t) = y(bridgeElement-2+i,localt)
  end do

  y(:,localt-2) = y(:,localt-1)
  y(:,localt-1) = y(:,localt)


! Finally the Great Time Loop:

do t=3,(nTimeSteps+1)
    y(1,localt) = bL1*y(1,localt-1)+bL2*y(2,localt-1)+bL3*y(3,localt-1)+bL4*y(1,localt-2)+bLF*0
    y(2,localt) = a1*(y(4,localt-1)-y(2,localt-1)+2*y(1,localt-1)) + a2*(y(3,localt-1)+y(1,localt-1)) &
 + a3*y(2,localt-1) + a4*y(2,localt-2) + a5*(y(3,localt-2)+y(1,localt-2))
    y(N-1,localt) = a1*(2*y(N,localt-1)-y(N-1,localt-1)+y(N-3,localt-1)) &
 + a2*(y(N,localt-1)+y(N-2,localt-1)) + a3*y(N-1,localt-1) + a4*y(N-1,localt-2) &
 + a5*(y(N,localt-2)+y(N-2,localt-2))
    y(N,localt) = bR1*y(N,localt-1) + br2*y(N-1,localt-1) + bR3*y(N-1,localt-1) + bR4*y(N,localt-2) &
 + bRF*0

    do i=3,(N-1)
      y(i,localt) = a1*(y(i+2,localt-1)+y(i-2,localt-1)) + a2*(y(i+1,localt-1) + y(i-1,localt-1)) &
 + a3*y(i,localt-1) + a4*y(i,localt-2) + a5*(y(i+1,localt-2)+y(i-1,localt-2)) &
 + deltaT*deltaT*N*stringForce(i)/Ms
    end do

    if (hammerDone.eqv..FALSE.) then
      call rigid_hammer(y(:,localt),hammerHeight,hammerMass,hammerVelocity,hammerForce, &
         deltaT,hammerDone,tension,bendingPrefactor, &
         k,b,rho,deltaX,vy,stringForce)
    end if
    


  if (t==dampT) then
    !!! Damped constants
    D = 1+b1Damped*deltaT
    nu = (2*b2Damped*deltaT)/(deltaX**2)
      
    a1 = (-(r**2)*mu)/D
    a2 = (r**2+4*(r**2)*mu+nu)/D
    a3 = (2-2*(r**2)-6*(r**2)*mu-2*nu)/D
    a4 = (-1+b1Damped*deltaT+2*nu)/D
    a5 = -nu/D
  
    b_R_denom = 1 + b1Damped*deltaT+zetaB*r
    bR1 = (2-2*r**2*mu - 2*r**2)/b_R_denom
    bR2 = (4*r**2*mu - 2*r**2)/b_R_denom
    bR3 = (-2*r**2*mu)/b_R_denom
    bR4 = (-1+b1Damped*deltaT + zetaB*r)/b_R_denom
    bRF = (deltaT**2/rho)/b_R_denom
    
    b_L_denom = 1 + b1Damped*deltaT+zetaL*r
    bL1 = (2-2*r**2*mu - 2*r**2)/b_L_denom
    bL2 = (4*r**2*mu - 2*r**2)/b_L_denom
    bL3 = (-2*r**2*mu)/b_L_denom
    bL4 = (-1+b1Damped*deltaT + zetaL*r)/b_L_denom
    bLF = (deltaT**2/rho)/b_L_denom
  end if

! Transfer the bridge element string amplitudes to bridgeY:
  do i=1,3
    bridgeY(i,t) = y(bridgeElement-2+i,localt)
  end do

  y(:,localt-2) = y(:,localt-1)
  y(:,localt-1) = y(:,localt)

end do

end subroutine


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine rigid_hammer(y,hammerHeight,hammerMass,hammerVelocity,hammerForce,deltaT, &
       hammerDone,tension,bendingPrefactor,k,b,rho,deltaX,vy,stringForce)

    real*8, intent(inout) :: y(:), vy(:), hammerHeight(:), hammerVelocity(:)
    real*8, intent(inout) :: hammerForce(:),stringForce(:)
    real*8, intent(in)    :: hammerMass, bendingPrefactor, rho, deltaT, deltaX, tension, k, b
    real*8                :: lastCompression(size(hammerHeight,1)), compression(size(hammerHeight,1)), &
         oldStringForce(size(y,1))
    logical, intent(out)  :: hammerDone
    logical               :: contact(size(y,1))
    integer               :: i, j, nContactElements

    nContactElements = size(y,1)
    hammerDone = .FALSE.
    contact = .TRUE.
    
    lastCompression = hammerHeight - y
    stringforce = 0
    
    forall (i=1:size(y,1),lastCompression(i).le.0)
       contact(i) = .FALSE.
       lastCompression(i) = 0
       nContactElements = nContactElements - 1
    end forall
    if (nContactElements.le.1) then
       hammerDone = .TRUE.
       return
    end if
    lastCompression = lastCompression**b
    
    
    hammerHeight = hammerHeight + hammerVelocity*deltaT - k*nContactElements*lastCompression*deltaT**2/(2*hammerMass)


    
    forall (j=1:size(y,1),contact(i).eqv..TRUE.)
       stringforce(i) = -(y(j-2)-4*y(j-1)-2*y(j)-4*y(j+1)+y(j+2))*(1+(y(j-1)-y(j+1))**2)
       stringforce(j) = -stringforce(j)*bendingPrefactor
       stringforce(j) = stringforce(j) - tension*(y(j-1)-2*y(j)+y(j+1))
    end forall
    
    y = 2*y - vy*deltaT + (k*lastCompression + stringforce)/(2*rho*deltaX)

!CALCULATEFORCES

    oldstringforce = stringforce
    compression = hammerHeight - y

    forall (i=1:size(y,1),compression(i).le.0)
       contact(i) = .FALSE.
       compression(i) = 0
       nContactElements = nContactElements - 1
    end forall
    if (nContactElements.le.1) then
       hammerDone = .TRUE.
       return
    end if
    compression = compression**b

    forall (i=1:size(y,1),contact(i).eqv..TRUE.)
       stringforce(j) = -(y(j-2)-4*y(j-1)-2*y(j)-4*y(j+1)+y(j+2))*(1+(y(j-1)-y(j+1))**2)
       stringforce(j) = -stringforce(j)*bendingPrefactor
       stringforce(j) = stringforce(j) - tension*(y(j-1)-2*y(j)+y(j+1))
    end forall
    
    hammerVelocity = hammerVelocity + k*nContactElements*deltaT*(lastcompression+compression)/(2*hammerMass)
    vy = vy + (k*(lastcompression+compression) + oldstringforce)*deltaT/(2*rho*deltaX)
    stringforce = k*compression
  end subroutine rigid_hammer
  




 

















  subroutine rigid_hammer_strike(y,hammerHeight,hammerMass,hammerVelocity,hammerForce,lastHammerVelocity,compression,deltaT, &
       profile,force,hammerDone,tension,bendingPrefactor,k,b,rho,deltaX)

    real*8, intent(inout) :: y(:), hammerHeight, hammerVelocity, lastHammerVelocity, compression, force(:)
    real*8, intent(inout) :: hammerForce
    real*8, intent(in)    :: profile(:), hammerMass, bendingPrefactor, rho, deltaT, deltaX, tension, k, b
    real*8                :: surface, oldforce(size(force,1)), oldcompression, comppowbp1, oldcomppowbp1, netforce, &
         oldHammerHeight, pstring, pstringmin1
    logical, intent(out)  :: hammerDone
    logical               :: contact(size(y,1))
    integer               :: i, j, nContactElements

    pstring = 0
    pstringmin1 = 0
    nContactElements = 0
    hammerDone = .FALSE.
    oldHammerHeight = hammerHeight
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
       force(j) = -(y(j-2)-4*y(j-1)-2*y(j)-4*y(j+1)+y(j+2))*(1+(y(j-1)-y(j+1))**2)
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
    comppowbp1 = comppowbp1 + oldcomppowbp1*(b+1)/k
    compression = comppowbp1**(1/(b+1))
    
    hammerHeight = hammerHeight + oldcompression - compression
    hammerForce = netforce - k*comppowbp1/compression

!    print *, hammerHeight, hammerVelocity, hammerForce, compression
  end subroutine rigid_hammer_strike



end module
