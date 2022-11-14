"""
用bfgs手动实现胞优化
"""
from ase.io import read
from ase.io import write
from deepmd.calculator import DP
from ase.spacegroup import crystal
from ase.io.trajectory import Trajectory
from ase.units import kB
from ase.optimize import BFGS
from ase.visualize import view
import numpy as np
from ase.eos import EquationOfState
from ase.eos import calculate_eos

# 计算产生的文件存放的文件夹
at_name = "attachment"

# 用 ase.spacegroup 创建 cgN
cgN_init = crystal(symbols=["N"], spacegroup=199, basis=[(0.085,0.085,0.085)], cellpar=[3.765, 3.765, 3.765, 90, 90, 90] )
cell = cgN_init.get_cell() # 获得胞长 -- 注意不能用cgN_init.cell，这返回的是view，而get_cell()返回的是拷贝
# 备份结构
cgN_init_copy = cgN_init.copy()

# 读取势
calc = DP(model="graph.pb")
# 将结构和势绑定
cgN_init.calc = calc

# 结构优化 - eos方法
# 将胞长手动在0.95 - 1.05的范围内变化，对其中每个结构做bfgs，取能量最小的值
traj = Trajectory(f"{at_name}/3_.traj", mode ='w')  # 储存每一步产生的结构
for x in np.linspace(0.97, 1.03, 15):
    cgN_init.set_cell(cell*x, scale_atoms=True)
    bfgs = BFGS(atoms=cgN_init)
    bfgs.run(fmax=0.01)
    traj.write(cgN_init)

# 读取结构
configs = read(f"{at_name}/3_.traj@0:")
volumes = [config.get_volume() for config in configs]
energies = [config.get_potential_energy() for config in configs]
# 计算状态方程
eos = EquationOfState(volumes, energies)
eos.plot(filename=f"{at_name}/3_eos.png")

"""
# 直接计算eos
eos = calculate_eos(cgN_init, npoints=10, eps=0.01)
eos.plot(filename=f"{at_name}/3_eos.png")
"""

# 改变胞的大小
v, e, B = eos.fit() # 体积，能量，体积模量
side_diff = (v/cgN_init.get_volume()) ** (1/3)
cgN_init.cell *= side_diff

print("开始位置 >>>")
print(cgN_init_copy.cell)
#print(cgN_init_copy.positions)
print("结束位置 >>>")
print(cgN_init.cell)
#print(cgN_init.positions)

# 储存CONTCAR
write(f"{at_name}/3_CONTCAR", images=cgN_init, format="vasp", parallel=True) # append参数：追加输入