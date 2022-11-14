"""
strainfilter结合bfgs可以进行胞优化
"""

from ase.io import read
from ase.io import write
from deepmd.calculator import DP
from ase.spacegroup import crystal
from ase.units import kB
from ase.optimize import BFGS
import numpy as np
from ase.constraints import StrainFilter

# 计算产生的文件存放的文件夹
at_name = "attachment"

# 用 ase.spacegroup 创建 cgN
cgN_init = crystal(symbols=["N"], spacegroup=199, basis=[(0.085,0.085,0.085)], cellpar=[3.765, 3.765, 3.765, 90, 90, 90] )
# 备份结构
cgN_init_copy = cgN_init.copy()

# 读取势
calc = DP(model="graph.pb")
# 将结构和势绑定
cgN_init.calc = calc

# 结构优化
sf = StrainFilter(atoms=cgN_init, mask=[1,1,1,1,1,1]) # mask是六个方向的应力都优化
bfgs = BFGS(atoms=sf, logfile=f"{at_name}/4_.log")
bfgs.run(fmax=0.001)

print("开始位置 >>>")
print(cgN_init_copy.cell)
print(cgN_init_copy.positions)
print("结束位置 >>>")
print(cgN_init.cell)
print(cgN_init.positions)

# 储存CONTCAR
write(f"{at_name}/4_CONTCAR", images=cgN_init, format="vasp", parallel=True) # append参数：追加输入