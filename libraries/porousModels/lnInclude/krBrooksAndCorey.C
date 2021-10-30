/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
  \\    /   O peration     |
  \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
  \\/     M anipulation  |
  -------------------------------------------------------------------------------
  License
  This file is part of OpenFOAM.

  OpenFOAM is free software: you can redistribute it and/or modify it
  under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
  FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
  for more details.

  You should have received a copy of the GNU General Public License
  along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

  \*---------------------------------------------------------------------------*/

#include "krBrooksAndCorey.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
  namespace relativePermeabilityModels
  {
    defineTypeNameAndDebug(krBrooksAndCorey, 0);

    addToRunTimeSelectionTable
    (
     relativePermeabilityModel,
     krBrooksAndCorey,
     dictionary
     );
  }
}




// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::relativePermeabilityModels::krBrooksAndCorey::krBrooksAndCorey
(
 const word& name,
 const dictionary& relativePermeabilityProperties,
 const volScalarField& Sb,
 const volScalarField& Sa
 )
  :
  relativePermeabilityModel(name, relativePermeabilityProperties,Sb,Sa),  
  Smin_
  (
      IOobject
      (
          "kr."+Sb_.name()+".r",
          Sb_.time().timeName(),
          Sb_.db(),
          IOobject::READ_IF_PRESENT,
          IOobject::NO_WRITE
      ),
      Sb.mesh(),
      relativePermeabilityProperties.subDict(typeName + "Coeffs").lookupOrDefault("kr."+Sb_.name()+".r",dimensionedScalar("kr."+Sb_.name()+".r",dimless,0))	  
  ),
  Smax_
  (
      IOobject
      (
          "kr."+Sa_.name()+".r",
          Sb_.time().timeName(),
          Sb_.db(),
          IOobject::READ_IF_PRESENT,
          IOobject::NO_WRITE
      ),
      Sb.mesh(),
      dimensionedScalar(Sa_.name()+"_r",dimless,1) - relativePermeabilityProperties.subDict(typeName + "Coeffs").lookupOrDefault("kr."+Sa_.name()+".r",dimensionedScalar("kr."+Sa_.name()+".r",dimless,0))  
  ),
  krBrooksAndCoreyCoeffs_(relativePermeabilityProperties.subDict(typeName + "Coeffs")),  
  nb_
  (
      IOobject
      (
          "kr.n."+Sb_.name(),
          Sb_.time().timeName(),
          Sb_.db(),
          IOobject::READ_IF_PRESENT,
          IOobject::NO_WRITE
      ),
      Sb.mesh(),
      krBrooksAndCoreyCoeffs_.lookupOrDefault<scalar>("kr.n."+Sb_.name(),0)
  ),
  na_
  (
      IOobject
      (
          "kr.n."+Sa_.name(),
          Sb_.time().timeName(),
          Sb_.db(),
          IOobject::READ_IF_PRESENT,
          IOobject::NO_WRITE
      ),
      Sb.mesh(),
      krBrooksAndCoreyCoeffs_.lookupOrDefault<scalar>("kr.n."+Sa_.name(),0)
  ),
  Se_
  (
   IOobject
   (
    name,
    Sb_.time().timeName(),
    Sb_.db(),
    IOobject::NO_READ,
    IOobject::NO_WRITE
    ),       
   Sb_
   ),
  kra_
  (
   IOobject
   (
    name,
    Sb_.time().timeName(),
    Sb_.db(),
    IOobject::NO_READ,
    IOobject::NO_WRITE
    ),
   Sb.mesh(),
   dimensionSet(0,0,0,0,0)
   ),
  krb_
  (
   IOobject
   (
    name,
    Sb_.time().timeName(),
    Sb_.db(),
    IOobject::NO_READ,
    IOobject::NO_WRITE
    ),       
   Sb.mesh(),
   dimensionSet(0,0,0,0,0)
  ),
  dkradS_
  (
   IOobject
   (
    name,
    Sb_.time().timeName(),
    Sb_.db(),
    IOobject::NO_READ,
    IOobject::NO_WRITE
    ),
   Sb.mesh(),
   dimensionSet(0,0,0,0,0)
   ),
  dkrbdS_
  (
   IOobject
   (
    name,
    Sb_.time().timeName(),
    Sb_.db(),
    IOobject::NO_READ,
    IOobject::NO_WRITE
    ),       
   Sb.mesh(),
   dimensionSet(0,0,0,0,0)
  ),
  krbmax_
  (
      IOobject
      (
          "kr.A."+Sb_.name(),
          Sb_.time().timeName(),
          Sb_.db(),
          IOobject::READ_IF_PRESENT,
          IOobject::NO_WRITE
      ),
      Sb.mesh(),
      krBrooksAndCoreyCoeffs_.lookupOrDefault<scalar>("kr.A."+Sb_.name(),1)
  ),
  kramax_
  (
      IOobject
      (
          "kr.A."+Sa_.name(),
          Sb_.time().timeName(),
          Sb_.db(),
          IOobject::READ_IF_PRESENT,
          IOobject::NO_WRITE
      ),
      Sb.mesh(),
      krBrooksAndCoreyCoeffs_.lookupOrDefault<scalar>("kr.A."+Sa_.name(),1.0)
  )
{
	
	Info << "krbmax_" << endl << krbmax_ << endl;
	Info << "kramax_" << endl << kramax_ << endl;
	
	
    if (gMin(nb_) <= 0 || gMin(na_) <= 0)
    {		
      FatalErrorIn
        (
         "in krBrooksAndCorey.C"
         )
        << "Relative permeability coefficient n equal or less than 0" 
        << exit(FatalError);
    }

}

// ************************************************************************* //
