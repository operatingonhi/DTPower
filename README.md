# DTPower

DTPower is a crank-angle-resolved simulation of a single-cylinder, four-stroke
spark-ignition engine. It treats the cylinder contents as a lumped thermodynamic
state with prescribed intake and exhaust boundary conditions. It models valve
flow, gas chemistry with Cantera, wall heat transfer with the Woschni correlation,
and piston friction. The executable runs a motored cycle to establish a
reference pressure trace, then runs a fired cycle.

## Build from source

The build below is verified on Linux with GCC 14, CMake 4.2, Cantera 3.2.0a4,
and yaml-cpp 0.8.0. You need a C++17 compiler, CMake 3.18 or newer, a build tool
such as Make or Ninja, `pkg-config`, the Cantera C++ development libraries and
mechanism data, yaml-cpp development libraries, and BLAS/LAPACK.

Follow the installation instructions for Cantera's C++ interface on your
platform. For dependencies in a custom prefix, add the directories containing
`cantera.pc` and `yaml-cpp.pc` to `PKG_CONFIG_PATH`.
[Cantera installation](https://cantera.org/stable/install/index.html) ·
[Cantera C++ build guide](https://cantera.org/stable/userguide/compiling-cxx.html)

```sh
git clone https://github.com/operatingonhi/DTPower.git
cd DTPower
pkg-config --modversion cantera yaml-cpp
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --parallel
ctest --test-dir build --output-on-failure
```

Set `-DBUILD_TESTING=OFF` at configure time for an executable-only build.
You can stage an installation with `cmake --install build --prefix "$PWD/build/install"`.
It includes the executable, example YAML files, README, and MIT license.
The installed executable uses the Cantera and other shared libraries from your
build environment. For Windows cross-compilation with MinGW, supply target
dependency paths through the toolchain file or CMake cache using
`CANTERA_INCLUDE_DIRS`, `CANTERA_LIBRARIES`, `YAML_CPP_INCLUDE_DIRS`, and
`YAML_CPP_LIBRARIES`. The release checks below were run on Linux.

## Run the methane–air example

From the repository root, copy the example into a writable run directory and
launch the program from there.

```sh
mkdir -p build/methane-air
cp examples/methane-air/*_parameters.yaml build/methane-air/
(cd build/methane-air && ../engine_sim)
```

The executable reads `engine_parameters.yaml` and `simulation_parameters.yaml`
from the working directory. The example uses Cantera's bundled `gri30.yaml`
mechanism. If it cannot be found, set `CANTERA_DATA` to the directory containing
it. Each run overwrites the following CSV files in the working directory.

| File | Contents |
| --- | --- |
| `motored_results.csv` | Crank-angle trace of the motored cycle |
| `engine_results.csv` | Fired-cycle pressure, temperature, volume, mass, heat transfer, and species mass fractions |
| `engine_performance.csv` | Mean effective pressures, power, efficiency, and an energy ledger |

## Configure a case

- Engine parameters include geometry in meters, engine speed in RPM, pressures in pascals,
  temperatures in kelvin, intake composition as Cantera mole ratios, and spark
  and valve settings in crank-angle degrees. Intake TDC is at 0°, compression
  TDC is at 360°, and positive `timing_btdc` advances the spark before 360°.
  [Example engine parameters](examples/methane-air/engine_parameters.yaml)
- Simulation parameters include the crank-angle step and output interval, heat-transfer and friction
  coefficients, and cycle range. Set `simulation.equilibrium_expansion: false`
  to freeze the composition after equilibrium ignition. The current YAML schema
  requires the `tolerances` section, whose pressure and temperature bounds are
  parsed but currently unused by the solver.
  [Example simulation parameters](examples/methane-air/simulation_parameters.yaml)

The `src/` directory contains the executable and model implementations,
`include/` holds their headers, `examples/` contains runnable inputs, and
`tests/` contains CTest regression sources.

## Validation status

The equilibrium regression checks internal energy during ignition, entropy and
element conservation during adiabatic expansion, energy balances for wall heat
transfer and discharge, frozen-composition behavior, and cycle reset.
The methane–air example exercises the full executable and produces an energy
ledger for inspection. Whole-cycle conservation, step-size convergence,
stability across operating conditions, comparison with reference data, and
handling of invalid or extreme inputs still need dedicated validation.

## License

DTPower is distributed under the [MIT License](LICENSE).
