import numpy as np
import ase
from ase.calculators.emt import EMT
from ase.io import read
from ase.io import write
from ase.spacegroup import crystal
from deepmd.calculator import DP
from ase.phonons import Phonons
import matplotlib.pyplot as plt

# 计算产生的文件存放的文件夹
at_name = "attachment"

# 用 ase.spacegroup 创建 cgN
cgN_init = crystal(symbols=["N"], spacegroup=199, basis=[(0.085,0.085,0.085)], cellpar=[3.765, 3.765, 3.765, 90, 90, 90])

# 读取势
calc = DP(model="graph.pb")

# 用DP势绑定结构
# cgN_init.calc = calc

# 声子计算
ph = Phonons(cgN_init, calc, supercell=(3,3,3), delta=0.001)
ph.run()

# 读取力并组装动力学矩阵
ph.read(acoustic=True)
ph.clean()

# 获得K_path
path = cgN_init.cell.get_bravais_lattice().bandpath(npoints=300)  # path.path >>> "GXMGRX,MR"
# 获得声子带结构 
bs = ph.get_band_structure(path)

# 计算dos
dos = ph.get_dos(kpts=(30, 30, 30)).sample_grid(npts=300, width=1e-3)

# 绘图
# plt.rcParams['font.sans-serif'] = ['Arial'] # 字体
fig = plt.figure(1, figsize=(7, 4), dpi=80)
ax = fig.add_axes([.12, .07, .67, .85])

emax = dos.get_energies().max()+0.002  # 纵坐标能量的最大值
bs.plot(ax=ax, emin=0.0, emax=emax, color="b", spin=None)
"""
energies = bs.energies  # energies.shape: 1,100,24  - 24条带
xcoords, label_xcoords, orig_labels = bs.get_labels() # xcoords.shape: 100
"""

dosax = fig.add_axes([.8, .07, .17, .85])
dosax.fill_between(x=dos.get_weights(), y1=dos.get_energies(), y2=0, color='grey',
                   edgecolor='k', lw=1)

dosax.set_ylim(0, emax)
dosax.set_yticks([])
dosax.set_xticks([])
dosax.set_xlabel("DOS", fontsize=12)

fig.savefig(f'{at_name}/8_cgN_phonon.png')


