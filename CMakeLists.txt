cmake_minimum_required(VERSION 3.8)

project(hybridPorousInterFoam)

set(CMAKE_CXX_STANDARD 14)

# set(CMAKE_OSX_ARCHITECTURES x86_64)
# set(CMAKE_OSX_ARCHITECTURES arm64)

include_directories(.)
include_directories(lnInclude)
include_directories(/opt/homebrew/include)
include_directories($ENV{FOAM_SRC}/OpenFOAM/lnInclude)
include_directories($ENV{FOAM_SRC}/OSspecific/POSIX/lnInclude)
include_directories($ENV{FOAM_SRC}/transportModels/lnInclude)
include_directories($ENV{FOAM_SRC}/twoPhaseModels/twoPhaseMixture/lnInclude)
include_directories($ENV{FOAM_SRC}/twoPhaseModels/interfaceProperties/lnInclude)
include_directories($ENV{FOAM_SRC}/twoPhaseModels/incompressibleTwoPhaseMixture/lnInclude)
include_directories($ENV{FOAM_SRC}/twoPhaseModels/immiscibleIncompressibleTwoPhaseMixture/lnInclude)
include_directories($ENV{FOAM_SRC}/twoPhaseModels/twoPhaseChange/lnInclude)
include_directories($ENV{FOAM_SRC}/MomentumTransportModels/momentumTransportModels/lnInclude)
include_directories($ENV{FOAM_SRC}/MomentumTransportModels/incompressible/lnInclude)
include_directories($ENV{FOAM_SRC}/finiteVolume/lnInclude)
include_directories($ENV{FOAM_SRC}/dynamicFvMesh/lnInclude)
include_directories($ENV{FOAM_SRC}/meshTools/lnInclude)
include_directories($ENV{FOAM_SRC}/sampling/lnInclude)
include_directories($ENV{FOAM_SRC}/topoChangerFvMesh/lnInclude)
include_directories(${PROJECT_SOURCE_DIR}/libraries/porousImmiscibleIncompressibleTwoPhaseMixture/lnInclude)
include_directories(${PROJECT_SOURCE_DIR}/libraries/porousInterfaceProperties/lnInclude)
include_directories(${PROJECT_SOURCE_DIR}/libraries/porousModels/lnInclude)

link_directories($ENV{FOAM_LIBBIN})
link_directories($ENV{FOAM_USER_LIBBIN})

if (APPLE AND EXISTS /opt/homebrew/lib)
    link_directories(/opt/homebrew/lib)
endif ()


add_definitions(
        -D$ENV{WM_ARCH}
        -DWM_ARCH_OPTION=$ENV{WM_ARCH_OPTION}
        -DWM_$ENV{WM_PRECISION_OPTION}
        -DWM_LABEL_SIZE=$ENV{WM_LABEL_SIZE}
        -DNoRepository -m64 -fPIC
)

find_package(nlohmann_json)

add_executable(${PROJECT_NAME} ${PROJECT_NAME}.C)

set(CMAKE_EXE_LINKER_FLAGS -Wl,-execute,-undefined,dynamic_lookup,-w)

target_link_libraries(${PROJECT_NAME}
        OpenFOAM
        immiscibleIncompressibleTwoPhaseMixture
        twoPhaseChangeModels
        momentumTransportModels
        incompressibleMomentumTransportModels
        finiteVolume
        dynamicFvMesh
        topoChangerFvMesh
        fvModels
        fvConstraints
        meshTools
        porousModels
        porousInterfaceProperties
        porousImmiscibleIncompressibleTwoPhaseMixture
        dl
        m        
        nlohmann_json::nlohmann_json)
