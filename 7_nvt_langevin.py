import numpy as np
from ase.io import read
from ase.io import write
from ase.build import make_supercell
from ase.spacegroup import crystal
from deepmd.calculator import DP
from ase.units import fs
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution
from ase.md.langevin import Langevin
from ase.io.vasp import write_vasp_xdatcar
from ase.constraints import StrainFilter
from ase.optimize import BFGS
from ase.md import MDLogger

# 计算产生的文件存放的文件夹
at_name = "attachment"

# 用 ase.spacegroup 创建 cgN
cgN_init = crystal(symbols=["N"], spacegroup=199, basis=[(0.085,0.085,0.085)], cellpar=[3.765, 3.765, 3.765, 90, 90, 90] )

# 读取势
calc = DP(model="graph.pb")

# 用DP势取计算结构
cgN_init.calc = calc

# 结构优化 - 防止NVT过程应力过大
sf = StrainFilter(cgN_init)
bfgs = BFGS(sf)
bfgs.run(fmax=1e-3)

# 扩胞
cgN_init = make_supercell(prim=cgN_init, P=3*np.eye(3))
cgN_init.calc = calc

# 分子动力学的初始温度
MaxwellBoltzmannDistribution(atoms = cgN_init, temperature_K=200) # 初始温度比正式跑的温度稍高，容易快速达到平衡态

# langevin - NVT
langevin = Langevin(atoms=cgN_init, timestep=1*fs, temperature_K=100, friction=1e-2, trajectory=f"{at_name}/7_md.traj", loginterval=1)

# 输出设置
langevin.attach(MDLogger(dyn=langevin, atoms=cgN_init, logfile=f"{at_name}/7_md.log", header=True, stress=True, peratom=True, mode="a"),
                interval=10)

langevin.run(steps=500)

# 将轨迹文件保存为xdatcar
write_vasp_xdatcar(f"{at_name}/7_xdatcar", images=read(f"{at_name}/7_md.traj@0:"), label="xdatcar for langevin md of cgN at 100 K")