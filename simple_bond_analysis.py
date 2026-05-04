"""
简化版键长分析脚本
专门针对YBa2Cu3O7结构的键长分析

使用方法:
python simple_bond_analysis.py YBa2Cu3O7_Ovac1_relaxed_MACE.cif
"""

from ase.io import read
import numpy as np
import sys


def analyze_bonds(cif_file):
    """
    分析CIF文件中的键长
    
    Parameters:
    -----------
    cif_file : str
        CIF文件路径
    """
    
    atoms = read(cif_file)
    
    print("="*60)
    print(f"键长分析: {cif_file}")
    print("="*60)
    print(f"化学式: {atoms.get_chemical_formula()}")
    print(f"原子数: {len(atoms)}")
    
    cell_lengths = atoms.get_cell_lengths_and_angles()
    print(f"晶胞参数:")
    print(f"  a = {cell_lengths[0]:.4f} Å")
    print(f"  b = {cell_lengths[1]:.4f} Å") 
    print(f"  c = {cell_lengths[2]:.4f} Å")
    print(f"  α = {cell_lengths[3]:.2f}°")
    print(f"  β = {cell_lengths[4]:.2f}°")
    print(f"  γ = {cell_lengths[5]:.2f}°")
    
    symbols = atoms.get_chemical_symbols()
    positions = atoms.get_positions()
    
    print(f"\n原子信息:")
    for i, (symbol, pos) in enumerate(zip(symbols, positions)):
        print(f"  {i+1:2d}. {symbol:2s} ({pos[0]:8.4f}, {pos[1]:8.4f}, {pos[2]:8.4f})")
    
    distances = atoms.get_all_distances(mic=True)
    
    bond_ranges = {
        'Cu-O': (1.8, 2.5),   # 铜氧键
        'Ba-O': (2.5, 3.2),   # 钡氧键  
        'Y-O':  (2.2, 2.8),   # 钇氧键
        'Cu-Cu': (2.5, 4.0),  # 铜铜距离
        'Ba-Ba': (3.5, 5.0),  # 钡钡距离
        'Y-Y':  (3.0, 5.0),   # 钇钇距离
        'O-O':  (2.5, 4.0),   # 氧氧距离
    }
    
    all_bonds = []
    
    for i in range(len(atoms)):
        for j in range(i + 1, len(atoms)):
            distance = distances[i, j]
            elem1, elem2 = symbols[i], symbols[j]
            
            bond_type = '-'.join(sorted([elem1, elem2]))
            
            is_reasonable = False
            for key, (min_d, max_d) in bond_ranges.items():
                if bond_type == key and min_d <= distance <= max_d:
                    is_reasonable = True
                    break
            
            if not is_reasonable and 0.8 <= distance <= 4.0:
                is_reasonable = True
            
            if is_reasonable:
                all_bonds.append((i+1, j+1, elem1, elem2, distance, bond_type))
    
    bond_groups = {}
    for bond in all_bonds:
        bond_type = bond[5]
        if bond_type not in bond_groups:
            bond_groups[bond_type] = []
        bond_groups[bond_type].append(bond)
    
    print(f"\n键长分析结果:")
    print("-"*80)
    
    for bond_type in sorted(bond_groups.keys()):
        bonds = bond_groups[bond_type]
        distances = [bond[4] for bond in bonds]
        
        print(f"\n{bond_type} 键 ({len(bonds)} 个):")
        print(f"  范围: {min(distances):.4f} - {max(distances):.4f} Å")
        print(f"  平均: {np.mean(distances):.4f} ± {np.std(distances):.4f} Å")
        
        bonds_sorted = sorted(bonds, key=lambda x: x[4])
        print(f"  详细信息:")
        for bond in bonds_sorted[:min(5, len(bonds))]:
            i, j, elem1, elem2, dist, _ = bond
            print(f"    原子{i:2d}({elem1}) - 原子{j:2d}({elem2}): {dist:.4f} Å")
        
        if len(bonds) > 5:
            print(f"    ... 还有 {len(bonds)-5} 个键")
    
    output_file = cif_file.replace('.cif', '_bonds.txt')
    with open(output_file, 'w') as f:
        f.write(f"键长分析结果: {cif_file}\n")
        f.write(f"化学式: {atoms.get_chemical_formula()}\n")
        f.write(f"原子数: {len(atoms)}\n\n")
        
        f.write(f"{'原子1':<6} {'原子2':<6} {'元素1':<4} {'元素2':<4} {'键长(Å)':<10} {'键类型':<8}\n")
        f.write("-"*50 + "\n")
        
        all_bonds_sorted = sorted(all_bonds, key=lambda x: x[4])
        for bond in all_bonds_sorted:
            i, j, elem1, elem2, dist, bond_type = bond
            f.write(f"{i:<6} {j:<6} {elem1:<4} {elem2:<4} {dist:<10.4f} {bond_type:<8}\n")
    
    print(f"\n详细结果已保存到: {output_file}")
    print("="*60)


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("使用方法: python simple_bond_analysis.py <cif_file>")
        print("示例: python simple_bond_analysis.py YBa2Cu3O7_Ovac1_relaxed_MACE.cif")
        sys.exit(1)
    
    cif_file = sys.argv[1]
    analyze_bonds(cif_file)
