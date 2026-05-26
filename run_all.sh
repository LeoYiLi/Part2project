#!/bin/bash
# 4个Frenkel构型分配到4张卡并行跑

CUDA_VISIBLE_DEVICES=0 python run_md.py YBa2Cu3O7_553_Frenkel.cif 0 &
CUDA_VISIBLE_DEVICES=1 python run_md.py YBa2Cu3O7_553_Frenkel2.cif 0 &
CUDA_VISIBLE_DEVICES=2 python run_md.py YBa2Cu3O7_553_Frenkel3.cif 0 &
CUDA_VISIBLE_DEVICES=3 python run_md.py YBa2Cu3O7_553_Frenkel4.cif 0 &

echo "所有任务已启动，等待完成..."
wait
echo "全部完成！"
