#!/bin/bash
# Quick Setup Script for Renalprog
# Run this after cloning the repository

echo "=========================================="
echo "Renalprog Installation Script"
echo "=========================================="
echo ""

# Check if mamba is available
if ! command -v mamba &> /dev/null; then
    echo "ERROR: mamba not found. Please install mamba/conda first."
    echo "Visit: https://github.com/conda-forge/miniforge"
    exit 1
fi

echo "Step 1: Creating conda environment..."
mamba create -n renalprog "python==3.9" -y

echo ""
echo "Step 2: Activating environment..."
eval "$(conda shell.bash hook)"
mamba activate renalprog

echo ""
echo "Step 3: Installing uv package manager..."
pip install uv

echo ""
echo "Step 4: Installing renalprog package..."
uv pip install -e .

echo ""
echo "Step 5: Installing test dependencies..."
uv pip install pytest pytest-cov

echo ""
echo "=========================================="
echo "Installation Complete!"
echo "=========================================="
echo ""
echo "Next steps:"
echo "1. Activate environment: mamba activate renalprog"
echo "2. Run tests: pytest tests/ -v"
echo ""
echo "To verify installation:"
echo "  python -c \"import renalprog; print(renalprog.__version__)\""
echo ""

