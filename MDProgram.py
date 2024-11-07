import numpy as np

def generate_non_lattice_chain(length, self_avoiding=False):
    """生成一个在三维空间中的随机行走路径，可选择是否避免自交。"""
    # 初始化起始点
    current_position = np.array([0.0, 0.0, 0.0])
    steps_and_lengths = [(0, 0.0)]  # 第一步（起点）为0，长度为0
    visited = {tuple(current_position)}  # 如果需要避免自交，用于记录已访问的位置
    
    for step_num in range(1, length + 1):
        # 生成一个随机方向（单位向量）
        theta = np.random.uniform(0, 2 * np.pi)
        phi = np.random.uniform(0, np.pi)
        
        # 将球坐标转换为笛卡尔坐标
        x = np.sin(phi) * np.cos(theta)
        y = np.sin(phi) * np.sin(theta)
        z = np.cos(phi)
        
        # 步长为1
        step = np.array([x, y, z])
        
        new_position = current_position + step
        
        if self_avoiding:
            # 检查是否需要避免自交
            if tuple(new_position) in visited:
                # 如果选择避免自交且新位置已被访问，则重新生成链
                return generate_non_lattice_chain(length, self_avoiding)
            else:
                visited.add(tuple(new_position))
        
        # 计算新位置到原点的距离
        distance = np.linalg.norm(new_position)
        steps_and_lengths.append((step_num, distance))
        current_position = new_position
    
    return np.array(steps_and_lengths)

def save_steps_and_lengths_to_file(steps_and_lengths, filename):
    """将步数和每一步的位置长度保存到 .dat 文件中。"""
    np.savetxt(filename, steps_and_lengths, header='step distance', comments='')

# 参数
chain_length = 100

# 分别生成两种情况的链
non_self_avoiding_steps_and_lengths = generate_non_lattice_chain(chain_length, self_avoiding=False)
self_avoiding_steps_and_lengths = generate_non_lattice_chain(chain_length, self_avoiding=True)

# 保存到文件
save_steps_and_lengths_to_file(non_self_avoiding_steps_and_lengths, 'non.dat')
save_steps_and_lengths_to_file(self_avoiding_steps_and_lengths, 'self.dat')

print("步数和位置长度已生成并保存到 .dat 文件：'non.dat' 和 'self.dat'")
