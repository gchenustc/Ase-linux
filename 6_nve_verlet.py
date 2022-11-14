import numpy as np
from ase.io import read
from ase.io import write
from ase.build import make_supercell
from ase.spacegroup import crystal
from deepmd.calculator import DP
from ase.units import fs
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution
from ase.md.verlet import VelocityVerlet

# 计算产生的文件存放的文件夹
at_name = "attachment"

# 用 ase.spacegroup 创建 cgN
cgN_init = crystal(symbols=["N"], spacegroup=199, basis=[(0.085,0.085,0.085)], cellpar=[3.765, 3.765, 3.765, 90, 90, 90] )

# 扩胞
cgN_init = make_supercell(prim=cgN_init, P=3*np.eye(3))

# 读取势
calc = DP(model="graph.pb")

# 用DP势取计算结构
cgN_init.calc = calc

# 分子动力学的初始温度
MaxwellBoltzmannDistribution(atoms = cgN_init, temperature_K=300)

# Verlet适合长时间的NVE模拟
# 在典型的NVE模拟中，温度将保持大致恒定，但如果发生显著的结构变化，则可能导致温度变化。如果有外部做功，温度可能会显著升高。
# loginterval关键字来指定写入轨迹的频率。loginterval关键字将应用于轨迹和日志文件。
verlet = VelocityVerlet(atoms=cgN_init, timestep=1*fs, logfile=f"{at_name}/6_.log", trajectory=f"{at_name}/6_.traj", loginterval=1)

verlet.run(1000)
