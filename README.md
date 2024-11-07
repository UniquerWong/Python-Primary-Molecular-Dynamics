 本文章介绍了用Python语言编写一个简单的初级分子动力学（MD）模拟框架程序，适用于计算化学的初步理解。

本程序用1000个粒子模拟体系从初始有序态到随机态的动力学过程。采用恒温恒容系综，恒温采用速度直接标度法。模拟空间是正立方体空间，要应用截断半径概念和周期边界条件。粒子为中性粒子且彼此之间只有范德华作用力，初始态粒子整齐均匀分布。

一、代码概述
1.粒子初始化：在一个三维盒子中生成并初始化N个粒子的初始位置和速度。
2.周期性边界条件（PBC）：在粒子出界时将其位置调整回盒子内，模拟无边界效果。
3.力计算：使用Lennard-Jones势计算粒子间的范德华力。
4.时间积分（Velocity-Verlet算法）：通过时间步长模拟粒子在力场下的运动。
5.温度校准：基于期望的温度对粒子的速度进行标定，以确保系统在给定温度下运行。
6.数据输出：将初始和最终结构输出到.gro文件中，可导入至VMD等分子可视化软件中进行进一步分析。

二、代码具体内容
1.参数设置
import numpy as np

# 参数设置
N = 1000         # 总粒子数
L = 10           # 系统的边长（立方体）
mass = 1         # 粒子的质量
sigma = 1        # 办模方程的参数，代表粒子间的距离尺度
epsilon = 1      # 办模方程的参数，代表势能深度
rcut = 2.5      # 力的计算范围
k_B_T = 1       # 温度标定（与玻尔兹曼常数的乘积）
time_steps = 1000  # 模拟的时间步数
dt = 0.01       # 时间步长

2. 初始化粒子位置和速度 (initialize_particles)
# 初始化粒子位置和速度
def initialize_particles():
    positions = np.zeros((N, 3))
    num_per_side = int(np.cbrt(N))
    spacing = L / num_per_side
    index = 0
    for i in range(num_per_side):
        for j in range(num_per_side):
            for k in range(num_per_side):
                positions[index] = np.array([i, j, k]) * spacing
                index += 1
                if index >= N:
                    break
            if index >= N:
                break
        if index >= N:
            break
    velocities = np.sqrt(2 * k_B_T / mass) * (np.random.rand(N, 3) - 0.5)
    return positions, velocities
代码中将粒子排列在一个立方体网格中，每个方向上的粒子数为

    num_per_side = 10
从而确保了粒子在初始时刻的均匀分布。



初始速度通过高斯分布生成，并与期望温度相关。

3. 设置周期性边界条件 (apply_pbc)
# 周期性边界条件
def apply_pbc(positions):
    return positions - L * np.round(positions / L)
使用周期性边界条件处理粒子的运动，以模拟无边界情况。若粒子超出盒子的范围，将其“反射”到盒子中，模拟无限大空间的效果。


4. 分子间范德华力的计算(calculate_forces)
范德华力用于计算粒子间相互作用的力。公式如下：




其中，表示深度势阱、是作用范围参数、是粒子间距离。

# 计算范德华力（矢量化）
def calculate_forces(positions):
    forces = np.zeros((N, 3))
    for i in range(N):
        disp = apply_pbc(positions[i] - positions)
        r2 = np.sum(disp ** 2, axis=1)
        mask = (r2 < rcut ** 2) & (r2 > 0)
        r2 = r2[mask]
        disp = disp[mask]
        r6_inv = (sigma ** 2 / r2) ** 3
        force_magnitudes = 48 * epsilon * (r6_inv ** 2 - 0.5 * r6_inv) / r2
        forces[i] += np.sum(force_magnitudes[:, None] * disp, axis=0)
        forces[mask] -= force_magnitudes[:, None] * disp
    return forces

5. 速度Verlet算法 (velocity_verlet)
程序使用速度Verlet算法更新粒子的运动状态，步骤如下：


更新位置：

         

计算新力：
         

更新速度：

      

这种算法由于使用当前位置和新位置的力来更新速度，使其具有较高的数值稳定性和能量守恒性。

# 速度Verlet算法
def velocity_verlet(positions, velocities, forces):
    positions += velocities * dt + 0.5 * forces * dt ** 2 / mass
    positions = apply_pbc(positions)
    new_forces = calculate_forces(positions)
    velocities += 0.5 * (forces + new_forces) * dt / mass
    return positions, velocities, new_forces
6. 温度标定速度 (rescale_velocities)

速度根据温度重新标定，以保证模拟在指定温度

k_B_T = 1
下进行。


通过缩放因子

v_scale
调整速度




其中，是速度平方的平均值。

# 温度标定速度
def rescale_velocities(velocities, k_B_T):
    v_scale = np.sqrt(k_B_T / (0.5 * np.mean(np.sum(velocities ** 2, axis=1))))
    velocities *= v_scale
    return velocities
7. GRO文件输出 (write_gro_file)
将粒子的位置信息写入

    initial_structure.gro
    final_structure.gro
文件，以标准格式输出，可供VMD等分子可视化软件进行进一步分析。

8. 主程序运行
# 主程序
positions, velocities = initialize_particles()
forces = calculate_forces(positions)
write_gro_file("initial_structure.gro", positions) #初始构象

for step in range(time_steps):
    positions, velocities, forces = velocity_verlet(positions, velocities, forces)
    velocities = rescale_velocities(velocities, k_B_T)

write_gro_file("final_structure.gro", positions) #最终构象


此分子动力学（MD）模拟框架为研究粒子系统在不同物理参数条件下的运动行为和相互作用特性提供了一个基础工具。通过调整模拟参数，如温度、截断半径、势参数、粒子数量和盒子尺寸等，研究人员可以探索不同条件对系统热力学性质、相变过程、扩散特性等方面的影响。此外，此框架可以进一步用于分析不同势函数对粒子相互作用的影响，为理解分子动力学提供了一个可拓展的程序架构。

​
