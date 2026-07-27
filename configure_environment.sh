#!/bin/bash
# Setup script for this repository

PYTHON=python3.13
VENV=.venv
BASE_DIR=$(pwd)
LIBRARIES=./libraries
MAVLINK=$LIBRARIES/mavlink
EIGEN=$LIBRARIES/mavlink
# Save base folders of header files
MAVLINK_GENERATED_HEADERS=$MAVLINK/generated/include/mavlink/v2.0

# First, initialize all submodules
echo "=== SUBMODULES ==="

echo "Initializing submodules..."
git submodule update --init --recursive
echo "Done!"

echo

# Then setup python venv for mavlink
echo "=== MAVLINK ==="

if [ -d $VENV ]; then
    echo "Found python virtual environment! Skipping installation."
else
    echo "Creating python virtual environment..."
    $PYTHON -m venv $VENV &>/dev/null
    echo "Done! Now activating..."
    source .venv/bin/activate
    echo "Done!"
fi

echo "Installing python dependencies..."

# Install custom dependencies for the project
python3 -m pip install -r requirements.txt &>/dev/null

# Install mavlink dependencies
cd $MAVLINK
python3 -m pip install -r pymavlink/requirements.txt &>/dev/null

# Build C code from mavlink submodule
echo "Done! Building MAVLINK codegen..."
python3 -m pymavlink.tools.mavgen \
            --lang=C \
            --wire-protocol=2.0 \
            --output=generated/include/mavlink/v2.0 \
            message_definitions/v1.0/common.xml
            #&>/dev/null
cd $BASE_DIR
echo "Done!"

echo

# Then install all other python dependencies for this project
echo "=== SPACECRAFT SIMULATOR ==="

python3 -m pip install -r requirements.txt

echo "Installed all python dependencies for Spacecraft Simulator!"
echo
