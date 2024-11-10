import numpy as np

# 参数设置
N = 1000
L = 10
mass = 1
sigma = 1
epsilon = 1
rcut = 2.5
dt = 0.01
k_B_T = 1
time_steps = 1000

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

# 周期性边界条件
def apply_pbc(positions):
    return positions - L * np.round(positions / L)

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

# 速度Verlet算法
def velocity_verlet(positions, velocities, forces):
    positions += velocities * dt + 0.5 * forces * dt ** 2 / mass
    positions = apply_pbc(positions)
    new_forces = calculate_forces(positions)
    velocities += 0.5 * (forces + new_forces) * dt / mass
    return positions, velocities, new_forces

# 温度标定速度
def rescale_velocities(velocities, k_B_T):
    v_scale = np.sqrt(k_B_T / (0.5 * np.mean(np.sum(velocities ** 2, axis=1))))
    velocities *= v_scale
    return velocities

# 输出gro文件
def write_gro_file(filename, positions):
    with open(filename, 'w') as f:
        f.write("MD Simulation Output\n")
        f.write(f"{N}\n")
        for i in range(N):
            particle_type = 'A' if i < N // 2 else 'B'
            f.write(f"{i + 1:5d}{particle_type:5s}{particle_type:5s}{i + 1:5d}{positions[i, 0]:8.3f}{positions[i, 1]:8.3f}{positions[i, 2]:8.3f}\n")
        f.write(f"   {L:10.5f}{L:10.5f}{L:10.5f}\n")

# 主程序
positions, velocities = initialize_particles()
forces = calculate_forces(positions)
write_gro_file("initial_structure.gro", positions)

for step in range(time_steps):
    positions, velocities, forces = velocity_verlet(positions, velocities, forces)
    velocities = rescale_velocities(velocities, k_B_T)

write_gro_file("final_structure.gro", positions)

