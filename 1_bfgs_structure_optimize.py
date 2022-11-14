"""
bfgs 只优化原子位置，不优化胞
"""
from ase.optimize import BFGS
from ase.io import read
from ase.io import write
from deepmd.calculator import DP

# 计算产生的文件存放的文件夹
at_name = "attachment"
# 读取结构
cgN_init = read("./POSCAR")
# 读取势
calc = DP(model="graph.pb")
# 将结构和势绑定
cgN_init.calc = calc

# 备份结构
cgN_init_copy = cgN_init.copy()

# bfgs
bfgs = BFGS(atoms=cgN_init, trajectory=f"{at_name}/1_bfgs.traj", restart=f"{at_name}/1_bfgs.pckl", logfile=f"{at_name}/1_bfgs.log")
bfgs.run(fmax=1E-2) # fmax 是力收敛的精度

print("开始位置 >>>")
print(cgN_init_copy.cell)
print(cgN_init_copy.positions)
print("结束位置 >>>")
print(cgN_init.cell)
print(cgN_init.positions)

# 储存CONTCAR
write(f"{at_name}/1_CONTCAR", images=cgN_init, format="vasp", parallel=True) # append参数：追加输入

