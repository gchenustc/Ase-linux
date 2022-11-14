"""
"""

from ase.io import read
from ase.io import write
from deepmd.calculator import DP
from ase.spacegroup import crystal
from ase.units import kB,GPa
from ase.optimize import BFGS
import numpy as np
from ase.constraints import ExpCellFilter
from ase.io.trajectory import Trajectory

# 计算产生的文件存放的文件夹
at_name = "attachment"

# 用 ase.spacegroup 创建 cgN
cgN_init = crystal(symbols=["N"], spacegroup=199, basis=[(0.085,0.085,0.085)], cellpar=[3.765, 3.765, 3.765, 90, 90, 90] )

# 读取势
calc = DP(model="graph.pb")

# 新建轨迹
traj = Trajectory(filename=f"{at_name}/5_.traj", mode="w")
# 结构优化
for press in [0,10,20,30,40,50]:
    cgN_init_copy = cgN_init.copy()
    # 将结构和势绑定
    cgN_init_copy.calc = calc
    
    # 施加压力
    ecf = ExpCellFilter(atoms=cgN_init_copy, mask=[1,1,1,1,1,1], scalar_pressure=press * GPa)  # 参数 hydrostatic_strain？
    # bfgs优化
    bfgs = BFGS(atoms=ecf, logfile=f"{at_name}/5_{press}GPa.log")
    bfgs.run(fmax=0.01)

    # 写入轨迹
    traj.write(cgN_init_copy)
    # 储存CONTCAR
    write(f"{at_name}/5_CONTCAR_{press}GPa", images=cgN_init_copy, format="vasp", parallel=True) # append参数：追加输入

