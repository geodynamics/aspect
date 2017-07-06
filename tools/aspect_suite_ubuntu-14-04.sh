#! /bin/bash

# 2015-06-18
# For Ubuntu 14.04 LTS
#
# This script creates the directory ${INSTALL_DIR} and uses it to store
# and build the dependencies for aspect.
# To speedup the build times, you can change ${MAKE_THREADS} to -jn
# where n is the number of threads your computer can run simultaneously.
# Beware that increasing the number of threads increases the memory requirement.
#
# cmake options taken mostly from aspect manual and dealii documentation.

INITIAL_CWD=$(pwd)

MAKE_THREADS="-j1"

INSTALL_DIR="aspect_suite"
ARCHIVE_DIR="${INSTALL_DIR}/archives"
BUILD_DIR="${INSTALL_DIR}/builds"
LIBRARY_DIR="${INSTALL_DIR}/libraries"
BINARY_DIR="${INSTALL_DIR}/binaries"

P4EST_TAR="p4est-1.1.tar.gz"
P4EST_DIR=$(basename --suffix=.tar.gz ${P4EST_TAR})
P4EST_URL="https://p4est.github.io/release/p4est-1.1.tar.gz"
P4EST_SHA1="ed8737d82ef4c97b9dfa2fd6e5134226f24c9b0b"
P4EST_SETUP="p4est-setup.sh"
P4EST_SETUP_URL="https://www.dealii.org/developer/external-libs/p4est-setup.sh"
P4EST_SETUP_SHA1="86ac6b7895c65ede2450b0ddbb51cc5a7d000d25"

TRILINOS_TAR="trilinos-release-12-0-1.tar.gz"
# trilinos has repetitive directory names
TRILINOS_DIR="Trilinos-"$(basename --suffix=.tar.gz ${TRILINOS_TAR})
TRILINOS_URL="https://github.com/trilinos/Trilinos/archive/trilinos-release-12-0-1.tar.gz"
TRILINOS_SHA1="e8c841aa0be9b17789b00c23530f2cd3a114091c"

DEALII_TAR="dealii-8.2.1.tar.gz"
DEALII_DIR=$(basename --suffix=.tar.gz ${DEALII_TAR})
DEALII_URL="https://github.com/dealii/dealii/releases/download/v8.2.1/dealii-8.2.1.tar.gz"
DEALII_SHA1="18a83feb7b2d9bb7c7b3d7721176a90aa505b1eb"

ASPECT_TAR="aspect-1.3.tar.gz"
ASPECT_DIR=$(basename --suffix=.tar.gz ${ASPECT_TAR})
ASPECT_URL="https://github.com/geodynamics/aspect/archive/v1.3.tar.gz"
ASPECT_SHA1="2c9fc07afdecc39a802a74b0516df90b2e786b07"


check_and_wget () {
	FILEPATH=$1
	FILESHA1=$2
	FILEURL=$3

	# Do we have the correct file?
	([[ -f "${FILEPATH}" ]] \
		&& (echo "${FILESHA1}" "${FILEPATH}" \
			| sha1sum --status --strict --warn --check -- >/dev/null)) \
	||	# get file
		((echo -e "\033[0;36m%%% Getting ${FILEURL}\033[0m" \
			&& wget --continue --output-document="${FILEPATH}" "${FILEURL}") \
		||	# download failed
			(echo -e "\033[0;31m%%% Failed to download $(basename ${FILEPATH})\033[0m" \
			&& exit 1)) \
		&&	# check file
			((echo "${FILESHA1}" "${FILEPATH}" \
				| sha1sum --status --strict --warn --check -- >/dev/null) \
		||	# check failed
			(echo -e "\033[0;31m%%% Mismatch of SHA1 for $(basename ${FILEPATH})\033[0m" \
			&& exit 1))
}

pretty_extract () {
	FILEPATH=$1

	# Has the archive already been extracted?
	[[ -f "${FILEPATH}.done" ]] \
	||	# extract archive
		(echo -e "\033[0;36m%%% Extracting $(basename ${FILEPATH})\033[0m" \
		&& tar --extract \
			--directory $(dirname "${FILEPATH}") \
			--file "${FILEPATH}" \
		&& touch "${FILEPATH}.done") \
	|| exit 1
}


make_needed_directories () {
	mkdir --verbose --parents \
		"${INSTALL_DIR}" \
		"${ARCHIVE_DIR}" \
		"${LIBRARY_DIR}" \
		"${BINARY_DIR}" \
		"${BUILD_DIR}" \
		"${BUILD_DIR}/${TRILINOS_DIR}" \
		"${BUILD_DIR}/${DEALII_DIR}" \
		"${BUILD_DIR}/${ASPECT_DIR}" \
	|| exit 1
}

download_pkgs() {
	# Do we have the right ubuntu pkgs?
	dpkg -s >/dev/null 2>&1 \
		build-essential \
		cmake \
		g++ \
		gcc \
		gfortran \
		git \
		libblas-dev \
		libhdf5-openmpi-dev \
		liblapack-dev \
		libopenmpi-dev \
		numdiff \
		openmpi-bin \
		zlib1g-dev \
	||	# get ubuntu pkgs
		(echo -e "\033[0;36m%%% Getting Ubuntu packages\033[0m" \
		&& echo -e "\n\033[1;31m%%% Requires sudo apt-get!\033[0m\n" \
		&& sudo apt-get install \
			build-essential \
			cmake \
			g++ \
			gcc \
			gfortran \
			git \
			libblas-dev \
			libhdf5-openmpi-dev \
			liblapack-dev \
			libopenmpi-dev \
			numdiff \
			openmpi-bin \
			zlib1g-dev \
		&& sudo -k \
		&& echo -e "\n\033[1;31m%%% sudo relinquished!\033[0m\n") \
	|| exit 1
}

build_p4est () {
	# Is p4est already built?
	[[ -f "${BUILD_DIR}/${P4EST_DIR}.done" ]] \
	&& return 0

	download_pkgs || exit 1

	make_needed_directories || exit 1

	# get p4est setup
	check_and_wget \
		"${ARCHIVE_DIR}/${P4EST_SETUP}" \
		"${P4EST_SETUP_SHA1}" \
		"${P4EST_SETUP_URL}" \
	|| exit 1

	# get p4est
	check_and_wget \
		"${ARCHIVE_DIR}/${P4EST_TAR}" \
		"${P4EST_SHA1}" \
		"${P4EST_URL}" \
	|| exit 1

	# extract p4est
	pretty_extract "${ARCHIVE_DIR}/${P4EST_TAR}" || exit 1

	# build p4est
	# cd # script creates a build folder in cwd
	(cd "${BUILD_DIR}" \
		&& echo -e "\033[0;36m%%% Running ${P4EST_SETUP}\033[0m" \
		&& bash \
			"${INITIAL_CWD}/${ARCHIVE_DIR}/${P4EST_SETUP}" \
			"${INITIAL_CWD}/${ARCHIVE_DIR}/${P4EST_DIR}" \
			"${INITIAL_CWD}/${LIBRARY_DIR}/${P4EST_DIR}" \
		&& cd "${INITIAL_CWD}" \
		&&	# clean p4est files
			# no p4est make clean
			rm -r \
				"${ARCHIVE_DIR}/${P4EST_DIR}" \
				"${ARCHIVE_DIR}/${P4EST_TAR}.done" \
		&& touch "${BUILD_DIR}/${P4EST_DIR}.done") \
	|| exit 1
}

build_trilinos () {
	# Is trilinos already built?
	[[ -f "${BUILD_DIR}/${TRILINOS_DIR}.done" ]] \
	&& return 0

	download_pkgs || exit 1

	make_needed_directories || exit 1

	# get trilinos
	check_and_wget \
		"${ARCHIVE_DIR}/${TRILINOS_TAR}" \
		"${TRILINOS_SHA1}" \
		"${TRILINOS_URL}" \
	|| exit 1

	# extract trilinos
	pretty_extract "${ARCHIVE_DIR}/${TRILINOS_TAR}" || exit 1

	# cd # cmake creates files in cwd
	(cd "${BUILD_DIR}/${TRILINOS_DIR}" \
		&& echo -e "\033[0;36m%%% Configuring ${TRILINOS_DIR}\033[0m" \
		&&	# configure trilinos
			# TPL_MPI_LIBRARIES # Trilinos cannot find openmpi libraries with Ubuntu openmpi 1.6.5-8
			cmake \
			-D Trilinos_ENABLE_Sacado=ON \
			-D Trilinos_ENABLE_Stratimikos=ON \
			-D CMAKE_BUILD_TYPE=RELEASE \
			-D CMAKE_CXX_FLAGS="-g -O3" \
			-D CMAKE_C_FLAGS="-g -O3" \
			-D CMAKE_FORTRAN_FLAGS="-g -O5" \
			-D Trilinos_EXTRA_LINK_FLAGS="-lgfortran" \
			-D CMAKE_VERBOSE_MAKEFILE=FALSE \
			-D Trilinos_VERBOSE_CONFIGURE=FALSE \
			-D TPL_ENABLE_MPI=ON \
			-D BUILD_SHARED_LIBS=ON \
			-D CMAKE_INSTALL_PREFIX="${INITIAL_CWD}/${LIBRARY_DIR}/${TRILINOS_DIR}" \
			-D TPL_MPI_LIBRARIES=/usr/lib/openmpi/include/mpi.h \
			"${INITIAL_CWD}/${ARCHIVE_DIR}/${TRILINOS_DIR}" \
		&&	# make trilinos
			echo -e "\033[0;36m%%% Making ${TRILINOS_DIR}\033[0m" \
		&& make ${MAKE_THREADS} install \
		&&	# clean trilinos files
			cd "${INITIAL_CWD}" \
		&& make -C "${BUILD_DIR}/${TRILINOS_DIR}" clean \
		&& rm -r \
			"${ARCHIVE_DIR}/${TRILINOS_DIR}" \
			"${ARCHIVE_DIR}/${TRILINOS_TAR}.done" \
		&& touch "${BUILD_DIR}/${TRILINOS_DIR}.done") \
	|| exit 1
}

build_dealii () {
	# Is dealii already built?
	[[ -f "${BUILD_DIR}/${DEALII_DIR}.done" ]] \
	&& return 0

	build_p4est || exit 1

	build_trilinos || exit 1

	# extract dealii
	pretty_extract "${ARCHIVE_DIR}/${DEALII_TAR}" || exit 1

	# cd # cmake creates files in cwd
	(cd "${BUILD_DIR}/${DEALII_DIR}" \
		&& echo -e "\033[0;36m%%% Configuring ${DEALII_DIR}\033[0m" \
		&&	# configure dealii
			# DEAL_II_LINKER_FLAGS # Linker cannot find openmpi libraries with Ubuntu openmpi 1.6.5-8
			cmake \
			-D DEAL_II_WITH_MPI=ON \
			-D DEAL_II_WITH_HDF5=ON \
			-D CMAKE_INSTALL_PREFIX="${INITIAL_CWD}/${LIBRARY_DIR}/${DEALII_DIR}" \
			-D TRILINOS_DIR="${INITIAL_CWD}/${LIBRARY_DIR}/${TRILINOS_DIR}" \
			-D P4EST_DIR="${INITIAL_CWD}/${LIBRARY_DIR}/${P4EST_DIR}" \
			-D CMAKE_VERBOSE_MAKEFILE=FALSE \
			-D DEAL_II_LINKER_FLAGS="-I /usr/lib/openmpi/include" \
			"${INITIAL_CWD}/${ARCHIVE_DIR}/${DEALII_DIR}" \
		&&	# make dealii
			echo -e "\033[0;36m%%% Making ${DEALII_DIR}\033[0m" \
		&& make ${MAKE_THREADS} install \
		&&	# clean dealii files
			cd "${INITIAL_CWD}" \
		&& make -C "${BUILD_DIR}/${DEALII_DIR}" clean \
		&& rm -r \
			"${ARCHIVE_DIR}/${DEALII_DIR}" \
			"${ARCHIVE_DIR}/${DEALII_TAR}.done" \
		&& touch "${BUILD_DIR}/${DEALII_DIR}.done") \
	|| exit 1
}

build_aspect () {
	# Is aspect already built?
	[[ -f "${BUILD_DIR}/${ASPECT_DIR}.done" ]] \
	&& return 0

	build_dealii || exit 1

	# extract aspect
	pretty_extract "${ARCHIVE_DIR}/${ASPECT_TAR}" || exit 1

	# cd # cmake creates files in cwd
	(cd "${BUILD_DIR}/${ASPECT_DIR}" \
		&&	# configure aspect
			echo -e "\033[0;36m%%% Configuring ${ASPECT_DIR}\033[0m" \
		&& cmake \
			-D DEAL_II_DIR="${INITIAL_CWD}/${LIBRARY_DIR}/${DEALII_DIR}" \
			-D EXECUTABLE_OUTPUT_PATH="${INITIAL_CWD}/${BINARY_DIR}" \
			-D CMAKE_BUILD_TYPE=RELEASE \
			"${INITIAL_CWD}/${ARCHIVE_DIR}/${ASPECT_DIR}" \
		&&	# make aspect
			echo -e "\033[0;36m%%% Making ${ASPECT_DIR}\033[0m" \
		&& make ${MAKE_THREADS} \
		&&	# clean aspect files
			# no make clean for aspect. it removes the binary
			cd "${INITIAL_CWD}" \
		&& rm -r \
			"${ARCHIVE_DIR}/${ASPECT_DIR}" \
			"${ARCHIVE_DIR}/${ASPECT_TAR}.done" \
		&& touch "${BUILD_DIR}/${ASPECT_DIR}.done") \
	|| exit 1
}


download_pkgs_main () {
	(download_pkgs \
		&& echo -e "\033[0;36m%%% Required Ubuntu packages installed.\033[0m") \
	|| exit 1
}

download_archives_main () {
	make_needed_directories || exit 1

	# get p4est setup
	check_and_wget \
		"${ARCHIVE_DIR}/${P4EST_SETUP}" \
		"${P4EST_SETUP_SHA1}" \
		"${P4EST_SETUP_URL}" \
	|| exit 1

	# get p4est
	check_and_wget \
		"${ARCHIVE_DIR}/${P4EST_TAR}" \
		"${P4EST_SHA1}" \
		"${P4EST_URL}" \
	|| exit 1

	# get trilinos
	check_and_wget \
		"${ARCHIVE_DIR}/${TRILINOS_TAR}" \
		"${TRILINOS_SHA1}" \
		"${TRILINOS_URL}" \
	|| exit 1

	# get dealii
	check_and_wget \
		"${ARCHIVE_DIR}/${DEALII_TAR}" \
		"${DEALII_SHA1}" \
		"${DEALII_URL}" \
	|| exit 1

	# get aspect
	check_and_wget \
		"${ARCHIVE_DIR}/${ASPECT_TAR}" \
		"${ASPECT_SHA1}" \
		"${ASPECT_URL}" \
	|| exit 1

	echo -e "\033[1;36m%%% All archives downloaded.\033[0m" || exit 1
}

download_all_main () {
	download_pkgs_main || exit 1

	download_archives_main || exit 1
}

build_p4est_main () {
	(build_p4est \
		&& echo -e "\033[1;36m%%% ${P4EST_DIR} built.\033[0m") \
	|| exit 1
}

build_trilinos_main () {
	(build_trilinos \
		&& echo -e "\033[1;36m%%% ${TRILINOS_DIR} built.\033[0m") \
	|| exit 1
}

build_dealii_main () {
	(build_dealii \
		&& echo -e "\033[1;36m%%% ${DEALII_DIR} built.\033[0m") \
	|| exit 1
}

build_aspect_main () {
	(build_aspect \
		&& echo -e "\033[1;36m%%% ${ASPECT_DIR} built.\033[0m") \
	|| exit 1
}

all_main () {
	(download_all_main \
		&& build_p4est_main \
		&& build_trilinos_main \
		&& build_dealii_main \
		&& build_aspect_main \
		&& echo -e "\n\n\033[1;35m%%% COMPLETE %%%\033[0m") \
	|| exit 1
}

help_main () {
	echo "This script downloads, configures, and compiles ASPECT and its"
	echo " required dependencies. It is written for Ubuntu 14.04 LTS."
	echo "NOTE: This script is based on snapshots of the software projects"
	echo " and does NOT get the latest versions."
	echo ""
	echo "Options:"
	echo "all - Performs every task required to build the ASPECT binary, which"
	echo "    will be placed in '${BINARY_DIR}'. This is a time and"
	echo "    memory intensive task."
	echo "download - Downloads both the required Ubuntu packages and the"
	echo "    tarballs (compressed archives) of each project."
	echo "download_packages - Makes an apt-get request to install the"
	echo "    necessary packages for building ASPECT and its dependencies."
	echo "    NOTE: This uses sudo."
	echo "download_archives - Downloads the tarballs of ASPECT and its"
	echo "    dependencies, checking the sha1s for accuracy. The downloads are"
	echo "    stored in '${ARCHIVE_DIR}'."
	echo "build_p4est - Builds p4est library files, placing them in"
	echo "    '${LIBRARY_DIR}/${P4EST_DIR}'."
	echo "build_trilinos - Builds Trilinos library files, placing them in"
	echo "    '${LIBRARY_DIR}/${TRILINOS_DIR}'. This is a"
	echo "    time and memory intensive task."
	echo "build_dealii - Builds deal.II library files, placing them in"
	echo "    '${LIBRARY_DIR}/${DEALII_DIR}'. If necessary, p4est and Trilinos"
	echo "    are also built. This is a time and memory intensive task."
	echo "build_aspect - Builds ASPECT binary, placing it in '${BINARY_DIR}'."
	echo "    If necessary, will also build all dependencies."
	echo ""
	echo "Projects:"
	echo "    ASPECT:   v1.3     https://aspect.dealii.org/"
	echo "    deal.II:  v8.2.1   https://dealii.org/"
	echo "    Trilinos: v12.0.1  https://trilinos.org/"
	echo "    p4est:    v1.1     http://p4est.org/"
}

menu_main () {
	OPTIONS="Options: all, download, download_packages, download_archives,\n build_p4est, build_trilinos, build_dealii, build_aspect, help (-h, --help)"

	if [[ -z $@ ]];
	then echo -e ${OPTIONS};
	elif ! [[ $@ == $1 ]];
	then echo "Error: Takes one input only" \
		&& echo -e ${OPTIONS} \
		&& exit 1;
	elif [[ $@ == "all" ]];
	then all_main || exit 1;
	elif [[ $@ == "download" ]];
	then download_all_main || exit 1;
	elif [[ $@ == "download_packages" ]];
	then download_pkgs_main || exit 1;
	elif [[ $@ == "download_archives" ]];
	then download_archives_main || exit 1;
	elif [[ $@ == "build_p4est" ]];
	then build_p4est_main || exit 1;
	elif [[ $@ == "build_trilinos" ]];
	then build_trilinos_main || exit 1;
	elif [[ $@ == "build_dealii" ]];
	then build_dealii_main || exit 1;
	elif [[ $@ == "build_aspect" ]];
	then build_aspect_main || exit 1;
	elif [[ $@ == "help" || $@ == "-h" || $@ == "--help" ]];
	then (echo -e ${OPTIONS} && help_main) || exit 1;
	else echo "Error: Not an accepted option" \
		&& echo -e ${OPTIONS} \
		&& exit 1;
	fi
}

menu_main $@ || exit 1


