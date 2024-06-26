set(CONFIG_ARGS "--host=cray")
set(HDF5_PATH "$ENV{HDF5_PATH}")
set(NETCDF_PATH "$ENV{NETCDF_PATH}")
set(PNETCDF_PATH "$ENV{PNETCDF_PATH}")
string(APPEND SLIBS " -L$ENV{NCAR_ROOT_HDF5}/lib -lhdf5_hl -lhdf5 -L$ENV{PNETCDF_PATH}/lib -lpnetcdf -L$ENV{NETCDF_PATH}/lib -lnetcdff -lnetcdf -L/glade/u/apps/derecho/23.09/spack/opt/spack/parallelio/2.6.2/cray-mpich/8.1.27/oneapi/2024.0.2/tjs5/lib -lpiof -lpioc -mkl ")
string(APPEND CFLAGS " -qopt-report -march=core-avx2 -mcmodel=large")
string(APPEND FFLAGS " -qopt-report -march=core-avx2 -mcmodel=large")
string(APPEND LDFLAGS "  -mcmodel=large -Wl,--copy-dt-needed-entries")
if (COMP_NAME STREQUAL gptl)
  string(APPEND CPPDEFS " -DHAVE_SLASHPROC")
endif()
if (COMP_NAME STREQUAL mpi-serial)
  string(APPEND CFLAGS " -std=c89 ")
endif()	
set(SCC icx)
set(SFC ifort)
