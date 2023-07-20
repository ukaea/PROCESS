#!/usr/bin/bash

echo "Converting Python PROCESS installation:"

PROCESS_PYTHON_INSTALL=$(python3 -c "import site, os; print(os.path.join(site.getsitepackages()[0], 'process'))" | awk '{$1=$1};1')
PYTHON_ABI_VERSION=$(python3 -c "import sysconfig; print(sysconfig.get_config_var('EXT_SUFFIX'))" | awk '{$1=$1};1')

echo "    PROCESS Install Folder:' ${PROCESS_PYTHON_INSTALL}"
echo "    Shared Object to be Converted: _fortran${PYTHON_ABI_VERSION}"

if [ ! -f "${PROCESS_PYTHON_INSTALL}/_fortran${PYTHON_ABI_VERSION}" ]
then
	echo "    ERROR: File '${PROCESS_PYTHON_INSTALL}/_fortran${PYTHON_ABI_VERSION}' not found, are you using the same Python version as the install?"
	exit 1
fi

for VAR in libprocess libgfortran
do
	if [ ! -f "${PROCESS_PYTHON_INSTALL}/lib/${VAR}.dylib" ]
	then
		echo "    ERROR: File '${PROCESS_PYTHON_INSTALL}/lib/${VAR}.dylib' could not be found."
		exit 1
	fi
	echo "    Running: install_name_tool -change @rpath/${VAR}.dylib ${PROCESS_PYTHON_INSTALL}/lib/${VAR}.dylib ${PROCESS_PYTHON_INSTALL}/_fortran${PYTHON_ABI_VERSION}"
	install_name_tool -change @rpath/${VAR}.dylib ${PROCESS_PYTHON_INSTALL}/lib/${VAR}.dylib ${PROCESS_PYTHON_INSTALL}/_fortran${PYTHON_ABI_VERSION}
done
