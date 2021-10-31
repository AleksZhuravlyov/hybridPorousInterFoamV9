/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
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

Application
    interFoam

Description
    Solver for 2 incompressible, isothermal immiscible fluids using a VOF
    (volume of fluid) phase-fraction based interface capturing approach,
    with optional mesh motion and mesh topology changes including adaptive
    re-meshing.

\*---------------------------------------------------------------------------*/

#include <set>
#include <map>
#include <cmath>
#include <string>
#include <iomanip>

#include <nlohmann/json.hpp>

using json = nlohmann::json;

#include "fvCFD.H"
#include "dynamicFvMesh.H"
#include "CMULES.H"
#include "EulerDdtScheme.H"
#include "localEulerDdtScheme.H"
#include "CrankNicolsonDdtScheme.H"
#include "subCycle.H"
#include "porousImmiscibleIncompressibleTwoPhaseMixture.H"
#include "noPhaseChange.H"
#include "noPhaseChange.H"
#include "kinematicMomentumTransportModel.H"
#include "pimpleControl.H"
#include "pressureReference.H"
#include "fvModels.H"
#include "fvConstraints.H"
#include "CorrectPhi.H"
#include "fvcSmooth.H"

#include "capillarityModel.H"
#include "relativePermeabilityModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[]) {
  #include "postProcess.H"

  #include "setRootCaseLists.H"
  #include "createTime.H"
  #include "createDynamicFvMesh.H"
  #include "initContinuityErrs.H"
  #include "createDyMControls.H"
  #include "createFields.H"
  #include "createFieldRefs.H"
  #include "createAlphaFluxes.H"
  #include "initCorrectPhi.H"
  #include "createUfIfPresent.H"

  turbulence->validate();

  if (!LTS) {
    #include "CourantNo.H"
    #include "setInitialDeltaT.H"
  }


  json interFoamData;

  std::map<std::string, double> timesUMgnAvs;
  std::map<std::string, double> timesAlphaAvs;
  std::map<std::string, double> timesAlphaUmgnAvs;
  std::map<std::string, double> timesFAvs;
  std::map<std::string, double> timesQInj;

  std::map<std::string, double> phasesNus;

  phasesNus["1"] = nu1.value();
  phasesNus["2"] = nu2.value();

  std::map<std::string, double> phasesRhos;
  phasesRhos["1"] = rho1.value();
  phasesRhos["2"] = rho2.value();

  IOdictionary transportPropertiesMy(IOobject("transportProperties",
                                              runTime.constant(),
                                              mesh,
                                              IOobject::MUST_READ_IF_MODIFIED,
                                              IOobject::NO_WRITE));

  dimensionedScalar sigma("sigma", dimensionSet(1, 0, -2, 0, 0, 0, 0), transportPropertiesMy.lookup("sigma"));

  const scalarField &volumeCell = mesh.V();

  auto volume = volumeCell * eps;

  auto totalVolume = gSum(volume);

  auto sigmaValue = sigma.value();
  if (Pstream::master()) {
    interFoamData["nus"] = phasesNus;
    interFoamData["rhos"] = phasesRhos;
    interFoamData["sigma"] = sigmaValue;
    interFoamData["V"] = totalVolume;
  }

  std::string inletBoundaryName = "bottom";

  const label &inletId = mesh.boundaryMesh().findPatchID(inletBoundaryName);

  // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
  Info << "\nStarting time loop\n" << endl;

  while (pimple.run(runTime)) {
    #include "readDyMControls.H"

    if (LTS) {
      #include "setRDeltaT.H"
    } else {
      #include "CourantNo.H"
      #include "alphaCourantNo.H"
      #include "setDeltaT.H"
    }

    runTime++;

    Info << "Time = " << runTime.timeName() << nl << endl;

    // --- Pressure-velocity PIMPLE corrector loop
    while (pimple.loop()) {
      if (pimple.firstPimpleIter() || moveMeshOuterCorrectors) {

        // Store divU from the previous mesh so that it can be mapped
        // and used in correctPhi to ensure the corrected phi has the
        // same divergence
        tmp<volScalarField> divU;

        if
            (
            correctPhi
            && !isType<twoPhaseChangeModels::noPhaseChange>(phaseChange)
            ) {
          // Construct and register divU for mapping
          divU = new volScalarField
              (
                  "divU0",
                  fvc::div(fvc::absolute(phi, U))
              );
        }

        fvModels.preUpdateMesh();


        mesh.update();

        if (mesh.changing()) {
          // Do not apply previous time-step mesh compression flux
          // if the mesh topology changed
          if (mesh.topoChanging()) {
            talphaPhi1Corr0.clear();
          }

          MRF.update();

          if (correctPhi) {
            // Calculate absolute flux
            // from the mapped surface velocity
            phi = mesh.Sf() & Uf();

            #include "correctPhi.H"

            // Make the flux relative to the mesh motion
            fvc::makeRelative(phi, U);
          }

          mixture.correct();

          #include "updateVariables.H"
          #include "updateDarcyVelocities.H"

          if (checkMeshCourantNo) {
            #include "meshCourantNo.H"
          }
        }
      }

      #include "alphaControls.H"
      #include "alphaEqnSubCycle.H"

      mixture.correct();

      #include "updateVariables.H"
      #include "UEqn.H"

      // --- Pressure corrector loop
      while (pimple.correct()) {
        #include "pEqn.H"
      }

      if (pimple.turbCorr()) {
        turbulence->correct();
      }
    }

    runTime.write();

    auto time = runTime.timeName();

    if (runTime.writeTime()) {
      timesAlphaAvs[time] = gAverage(alpha1);
      timesUMgnAvs[time] = gAverage(mag(U)->v());
      timesAlphaUmgnAvs[time] = gAverage((mag(U) * alpha1)->v());

      timesFAvs[time] = timesAlphaUmgnAvs[time] / timesUMgnAvs[time];

      timesQInj[time] = gSum(*(((U.boundaryField()[inletId]) & (mesh.Sf().boundaryField()[inletId])).ptr()));


      Info << "achievedAlpha " << timesAlphaAvs[runTime.timeName()] << endl;


      if (Pstream::master()) {
        interFoamData["times_u_mgn_avs"] = timesUMgnAvs;
        interFoamData["times_alpha_avs"] = timesAlphaAvs;
        interFoamData["times_alpha_u_mgn_avs"] = timesAlphaUmgnAvs;
        interFoamData["times_F_avs"] = timesFAvs;
        interFoamData["times_Q_inj"] = timesQInj;

        std::ofstream ofs("results_" + time + ".json");
        ofs << std::setw(4) << interFoamData << std::endl;
      }

    }


    Info << "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
         << "  ClockTime = " << runTime.elapsedClockTime() << " s"
         << nl << endl;
  }

  if (Pstream::master()) {
    interFoamData["times_u_mgn_avs"] = timesUMgnAvs;
    interFoamData["times_alpha_avs"] = timesAlphaAvs;
    interFoamData["times_alpha_u_mgn_avs"] = timesAlphaUmgnAvs;
    interFoamData["times_F_avs"] = timesFAvs;
    interFoamData["times_Q_inj"] = timesQInj;

    std::ofstream ofs("results.json");
    ofs << std::setw(4) << interFoamData << std::endl;
  }

  Info << "End\n" << endl;

  return 0;
}


// ************************************************************************* //
