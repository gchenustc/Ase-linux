"""
basinHopping: 只优化原子位置，不优化胞，可以提供温度克服势垒
"""
from ase.io import read
from ase.io import write
from deepmd.calculator import DP
from ase.spacegroup import crystal
from ase.build import make_supercell
from ase.units import kB
from ase.optimize import BFGS
from ase.optimize.basin import BasinHopping
from ase.visualize import view
import numpy as np

# 计算产生的文件存放的文件夹
at_name = "attachment"

# 用 ase.spacegroup 创建 cgN
cgN_init = crystal(symbols=["N"], spacegroup=199, basis=[(0.085,0.085,0.085)], cellpar=[3.765, 3.765, 3.765, 90, 90, 90] )
# 扩胞 
P = 3*np.eye(3)
cgN_init = make_supercell(prim=cgN_init, P=P, )
# 备份结构
cgN_init_copy = cgN_init.copy()

# 读取势
calc = DP(model="graph.pb")
# 将结构和势绑定
cgN_init.calc = calc

# 全局结构优化
BasinHopping(atoms=cgN_init, temperature=100*kB, optimizer=BFGS, dr=0.5, fmax=0.1) # dr:步长


print("开始位置 >>>")
print(cgN_init_copy.cell)
print(cgN_init_copy.positions)
print("结束位置 >>>")
print(cgN_init.cell)
print(cgN_init.positions)

# 储存CONTCAR
write(f"{at_name}/2_CONTCAR", images=cgN_init, format="vasp", parallel=True) # append参数：追加输入