export HDF5_INSTALL_PREFIX=$FIREMODELS/"$FDS_BUILD_TARGET"_libs/hdf5/$HDF5_VERSION

if [[ "$FDS_BUILD_TARGET" == "ompi_intel"* ]]; then
../configure --enable-fortran --enable-parallel CXX=$(which mpicxx) CC=$(which mpicc) FC=$(which mpifort) --prefix=$HDF5_INSTALL_PREFIX
elif [[ "$FDS_BUILD_TARGET" == "ompi_gnu"* ]]; then
../configure --enable-fortran --enable-parallel CXX="$(which mpicxx)" CC="$(which mpicc)" FC="$(which mpifort)" LDFLAGS="-L$(dirname $(which mpifort))/../lib" --prefix=$HDF5_INSTALL_PREFIX
elif [[ "$FDS_BUILD_TARGET" == "impi_intel"* ]]; then
../configure  --disable-shared --enable-fortran --enable-parallel CXX="$(which mpiicpc) -cc=$(which icpx)" CC="$(which mpiicc) -cc=$(which icx)" FC="$(which mpiifx) -fc=$(which ifx)" --prefix=$HDF5_INSTALL_PREFIX
else
../configure --enable-fortran --enable-parallel CXX="$(which mpicxx)" CC="$(which mpicc)" FC="$(which mpifort)" --prefix=$HDF5_INSTALL_PREFIX
fi

make install



