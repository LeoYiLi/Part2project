#!/bin/bash
# 一键部署 YBCO MD 环境（配合 PyTorch 镜像使用）
# 使用方法: bash setup.sh

set -e

echo "=== 安装 MACE 和依赖 ==="
pip install mace-torch ase scipy matplotlib spglib

echo "=== 验证安装 ==="
python -c "import torch; print(f'PyTorch: {torch.__version__}, CUDA: {torch.cuda.is_available()}')"
python -c "from mace.calculators import MACECalculator; print('MACE: OK')"

echo "=== 完成！==="
echo "直接运行: bash run_all.sh"
