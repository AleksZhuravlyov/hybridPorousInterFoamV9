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

#include "pcBrooksAndCorey.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
  namespace capillarityModels
  {
    defineTypeNameAndDebug(pcBrooksAndCorey, 0);

    addToRunTimeSelectionTable
    (
     capillarityModel,
     pcBrooksAndCorey,
     dictionary
     );
  }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::capillarityModels::pcBrooksAndCorey::pcBrooksAndCorey
(
 const word& name,
 const dictionary& capillarityProperties,
 const volScalarField& Sb
 )
    :
  capillarityModel(name, capillarityProperties,Sb),	
  pcBrooksAndCoreyCoeffs_(capillarityProperties.subDict(typeName + "Coeffs")),
  Sminpc_
  (
      IOobject
      (
          "pc."+Sb_.name()+".min",
          Sb_.time().timeName(),
          Sb_.db(),
          IOobject::READ_IF_PRESENT,
          IOobject::NO_WRITE
      ),
      Sb.mesh(),
      pcBrooksAndCoreyCoeffs_.lookupOrDefault("pc."+Sb_.name()+".min",capillarityProperties.lookupOrDefault("pc."+Sb_.name()+".min",dimensionedScalar("pc."+Sb_.name()+".min",dimless,0)))
  ),
  Smaxpc_
  (
      IOobject
      (
          "pc."+Sb_.name()+".max",
          Sb_.time().timeName(),
          Sb_.db(),
          IOobject::READ_IF_PRESENT,
          IOobject::NO_WRITE
      ),
      Sb.mesh(),
      pcBrooksAndCoreyCoeffs_.lookupOrDefault("pc."+Sb_.name()+".max",capillarityProperties.lookupOrDefault("pc."+Sb_.name()+".max",dimensionedScalar("pc."+Sb_.name()+".max",dimless,0)))
  ),
  pc0_
  (
      IOobject
      (
          "pc.A",
          Sb_.time().timeName(),
          Sb_.db(),
          IOobject::READ_IF_PRESENT,
          IOobject::NO_WRITE
      ),
      Sb.mesh(),
      pcBrooksAndCoreyCoeffs_.lookupOrDefault("pc.A",dimensionedScalar("pc.A",dimless,0))
  ),
  alpha_
  (
      IOobject
      (
          "pc.n",
          Sb_.time().timeName(),
          Sb_.db(),
          IOobject::READ_IF_PRESENT,
          IOobject::NO_WRITE
      ),
      Sb.mesh(),
      pcBrooksAndCoreyCoeffs_.lookupOrDefault<scalar>("pc.n",0)
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
  pc_
  (
   IOobject
   (
    name,
    Sb.time().timeName(),
    Sb.db(),
    IOobject::NO_READ,
    IOobject::NO_WRITE
    ),       
   Sb.mesh(),
   dimensionSet(1,-1,-2,0,0,0,0)
   ),
  dpcdS_
  (
   IOobject
   (
    name,
    Sb.time().timeName(),
    Sb.db(),
    IOobject::NO_READ,
    IOobject::NO_WRITE
    ),       
   Sb.mesh(),
   dimensionSet(1,-1,-2,0,0,0,0)
   )
{
	
    if (gMin(alpha_) == 0) FatalErrorIn("Foam::capillarityModels::pcBrooksAndCorey::pcBrooksAndCorey") << "alpha = 0 in pcBrooksAndCorey" << abort(FatalError);
}

// ************************************************************************* //
