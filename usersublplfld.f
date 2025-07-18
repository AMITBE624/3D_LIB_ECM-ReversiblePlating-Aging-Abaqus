c
c Compute lithium plating currents at microscale
c used in battery - ECHEM Porous Electrode Theory
c
      subroutine uechemecdlpl(
c     OUTPUT
     *     ecdlpl,         
     *     dEcdlpl_dCSurfNorm, dEcdlpl_dT, dEcdlpl_dSol, dEcdlpl_dPress,
     *     overpotEta_lpl,
c     INPUT
     *     filmResis, chargeDischarge,
     *     jmac,jMatyp,matlayo,laccFlag,
     *     cmname, cParticleName,
     *     kstep,lFlags, kinc, time, dtime, noel, npt, 
     *     ntens, ndi, nshr, props, nprops,
     *     predef, dpredef, npredef,
     *     statev, statevOld, nstatev,
     *     Temp, dTemp, Press, coords,
     *     solutionvar, dsolutionvar, nSolutionVar,
     *     cSurfNorm, cAvgNorm, ecdTotalOld, 
     *     SEILPLFilmThick,concLPL,
     *     i_array,niarray,r_array,nrarray,c_array,ncarray)
c
c default Abaqus include file to declare variables based on "implicit" fortran definitions.
      INCLUDE 'ABA_PARAM.INC'
c Below is to access property and parameter table.
      include 'aba_tcs_param.inc' 
c below is to get universal constants.
      include 'aba_phcon_param.inc'
c below is for location of various solution variables.
      include 'aba_echemsolutionvarloc.inc'
c
      parameter (one = 1.0d0, zero=0.0d0, tol = 1e-6)
c
      dimension i_array(*), r_array(*), jmatyp(*),jmac(*), lFlags(*)
      dimension time(2), props(nprops), predef(npredef), dpredef(npredef),
     *     statev(nstatev), statevold(nstatev), coords(*),
     *     solutionvar(nsolutionvar), dsolutionvar(nsolutionvar),
     *     dEcdlpl_dSol(*)
      character*80 cmname, c_array(NCARRAY), cParticleName
c
      parameter (maxParamsEChem = 30)

c     for parameter table
      character*80 cParams(maxParamsEChem)
      dimension iParamDataTypes(maxParamsEChem), 
     *    rParams(maxParamsEChem), iParams(maxParamsEChem)
      character*80 cParamsParticles(maxParamsEChem)
      character*80 tableName, c_regionname, cParticleNameUsub, c_particlenameRead
c 
      tempnew = temp
      phiS = solutionVar(i_sol_PhiS)
      phiE = solutionVar(i_sol_PhiE)
      aIonConc = solutionvar(i_sol_Ce)
c
cccc Get the physical constants      
      R_UnivConst = zero
      t_abs0 = zero
      FaradayConst  = zero
      jError = 0
      call getphysicalconstant(j_phcon_UnivGas, R_UnivConst, jError)
      jError = 0
      call getphysicalconstant(j_phcon_Abs0,  t_abs0, jError)
      jError = 0
      call getphysicalconstant(j_phcon_faraday, FaradayConst, jError)     
c
      oneOverT = zero
      if(abs(tempnew-t_abs0).gt.tol)oneOverT = one/(tempnew-t_abs0)
      oneOverRT = one/R_UnivConst * oneOverT 

c  Starting from electrode level to particle and layer level.
      tableName = 'ABQ_EChemPET_Electrode_Definition'
      jError = 0
      numParams = 0
      call getParameterTable(tableName, 
     *     numParams, iParamDataTypes, iParams, rParams, cParams, jError)
      c_regionname   = cParams(1)

c Check if particles are defined.
c    
      tableName = 'ABQ_EChemPET_'//TRIM(c_regionname)//'_Particles'
      jError = 0
      numRows=0
      numParticles = 0
      call queryParameterTable(tableName, numParams, numRows, jError)

      if (jError.eq.0) then 
         numParticles = numRows 
      end if
      
      numParamsParticles = 0
      if(numParticles.gt.0)then
         jerror = 0
         call getParameterTable(tableName, 
     *        numParamsParticles, iParamDataTypes, iParams,
     *        rParams, cParamsParticles, jError)
      endif

      do kParticle = 1, numParticles
         iParticleStart = (kParticle-1)*numParamsParticles
         c_particlenameRead  = cParamsParticles(iParticleStart+1)         

         if(c_particlenameRead .eq. cParticleName)then
            cParticleNameUsub = 'ABQ_EChemPET_'//TRIM(c_regionname)
     &           //'_' //TRIM(cParticleName)

c           Read properties for LPL current definitions -- START
            tableName = TRIM(cParticleNameUsub)//'_LPL'
            call getParameterTable(tableName, 
     *           numParams, iParamDataTypes, iParams, rParams, cParams, jError)

            alphacLpl = rParams(1)
            alphaaLpl = (1-alphacLpl)
            concLplInit =  rParams(2)
         endif
      end do

c  Compute the lithium plating current
c-----------------------------------------------------------------------------Plating/Stripping
c--      cfIrev: rate constant to define LPL irev. concentration  - eq(15) in paper <>
      cfIrev = 32E-7
c
      aI0_LPL = 0.076289949
c
c-- Reference of the following equation	  
      eq = log(aIonConc/conclpl)*((tempnew-t_abs0)*R_UnivConst)/FaradayConst
	  
      overpotEta_LPL = phiS - phiE - filmResis * ecdTotalOld - eq
      gammaLPL1 = alphaaLpl * FaradayConst * overpotEta_LPL * oneOverRT
      gammaLPL2 = - alphacLpl * FaradayConst * overpotEta_LPL * oneOverRT
	  
      stepN = (kstep - 2) * one
      rem=mod(stepN,4.)
      if(statevOld(1).le.zero)then
         arLi = concLplInit
         statev(1) = arLi
      else
         arLi = statevOld(1)
         rLi = max((conclpl - arLi),zero)

         swtch = 0
         if(rLi.gt.zero)then
            swtch = one
         endif	
         statev(1) = arLi + cfIrev * dtime * swtch * arLi
      endif

      if((rem.eq.zero).or.(rem.eq.3))then
         s2 = 0
         s1 = rLi/aIonConc
      else
         s2 = aIonConc/conclpl
         s1 = 0
      endif
      expgammaLPL1 = s1*exp(gammaLPL1)
      expgammaLPL2 = s2*exp(gammaLPL2)
      
      ecdlpl = aI0_LPL * (expgammaLPL1-expgammaLPL2)
      
      deqaIonConc = ((tempnew-t_abs0)*R_UnivConst)/(FaradayConst*aIonConc)
      deqT = log(aIonConc/conclpl)*R_UnivConst/FaradayConst
c
c  Linearizations for LPL current wrt different fields.
c
      gammaCoeff1 = alphaaLpl * FaradayConst * oneOverRT * overpotEta_LPL
      delgammaLPL_delPhiS1 = alphaaLpl * FaradayConst * oneOverRT
      delgammaLPL_delT1 = -gammaCoeff1 * oneOverT + 
     &     alphaaLpl * FaradayConst * oneOverRT * deqT
      delgammaLPL_delaIonConc1 = - alphaaLpl * FaradayConst * oneOverRT * deqaIonConc
	  
      gammaCoeff2 = -alphacLpl * FaradayConst * oneOverRT * overpotEta_LPL
      delgammaLPL_delPhiS2 = -alphacLpl * FaradayConst * oneOverRT
      delgammaLPL_delT2 = -gammaCoeff2 * oneOverT - 
     &     alphacLpl * FaradayConst * oneOverRT * deqT
      delgammaLPL_delaIonConc2 = alphacLpl * FaradayConst * oneOverRT * deqaIonConc
c   _delPhiS
      dEcdlpl_dSol(i_sol_phis) = aI0_LPL *(expgammaLPL1 * delgammaLPL_delPhiS1 -
     &                                     expgammaLPL2 * delgammaLPL_delPhiS2)
          
c   _delPhiE
      dEcdlpl_dSol(i_sol_phie) = - dEcdlpl_dSol(i_sol_phis)

c   _delCe
      T1 = aI0_LPL * (expgammaLPL1* delgammaLPL_delaIonConc1 - 
     &                expgammaLPL2 * delgammaLPL_delaIonConc2)
      T2 = aI0_LPL * (exp(gammaLPL1)*rLi/(aIonConc**2) + exp(gammaLPL2)/conclpl)
      dEcdlpl_dSol(i_sol_ce) = T1 - T2
      dEcdlpl_dT1 = ( expgammaLPL1 * delgammaLPL_delT1 *  aI0_LPL )
      dEcdlpl_dT2 = ( expgammaLPL2 * delgammaLPL_delT2 *  aI0_LPL )
      dEcdlpl_dT=dEcdlpl_dT1 - dEcdlpl_dT2	  

      return
      end
c -------------------------------
c
c Subroutine USDFLD
c 
c -------------------------------
      SUBROUTINE USDFLD(FIELD,STATEV,PNEWDT,DIRECT,T,CELENT,
     $     TIME,DTIME,CMNAME,ORNAME,NFIELD,NSTATV,NOEL,NPT,LAYER,
     $     KSPT,KSTEP,KINC,NDI,NSHR,COORD,JMAC,JMATYP,MATLAYO,LACCFLA)
C
      INCLUDE 'ABA_PARAM.INC'
C
      CHARACTER*80 CMNAME,ORNAME
      CHARACTER*3  FLGRAY(15)

      double precision ECDTA, Capacity, Area, Ampere, num_layer
      DIMENSION FIELD(NFIELD),STATEV(NSTATV),DIRECT(3,3),
     $     T(3,3),TIME(2)
      DIMENSION ARRAY(15),JARRAY(15),JMAC(*),JMATYP(*),COORD(*)

      CALL GETVRM('CONCSEI_1',ARRAY,JARRAY,FLGRAY,JRCD,JMAC,JMATYP,
     $     MATLAYO,LACCFLA)
   
      CONCSEI = ARRAY(1) 

      FIELD(1)  = CONCSEI + STATEV(1)
      STATEV(2) = CONCSEI + STATEV(1)
C
      RETURN
      END
