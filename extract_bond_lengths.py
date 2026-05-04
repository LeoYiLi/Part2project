"""
简单的键长提取脚本
从CIF文件中计算并输出所有原子间的键长

使用方法:
python extract_bond_lengths.py structure.cif
"""

from ase.io import read
import numpy as np
import sys
import argparse


def get_bond_lengths(atoms, max_distance=4.0, min_distance=0.5):
    """
    计算结构中所有原子间的键长
    
    Parameters:
    -----------
    atoms : ase.Atoms
        原子结构
    max_distance : float
        最大键长阈值 (Å)
    min_distance : float
        最小键长阈值 (Å)
    
    Returns:
    --------
    bonds : list
        键长信息列表 [(atom1_idx, atom2_idx, element1, element2, distance), ...]
    """
    
    bonds = []
    symbols = atoms.get_chemical_symbols()
    positions = atoms.get_positions()
    
    distances = atoms.get_all_distances(mic=True)  # mic=True 考虑周期性边界条件
    
    for i in range(len(atoms)):
        for j in range(i + 1, len(atoms)):
            distance = distances[i, j]
            
            if min_distance <= distance <= max_distance:
                bonds.append((
                    i, j,
                    symbols[i], symbols[j],
                    distance
                ))
    
    return bonds


def print_bond_lengths(bonds, sort_by='distance'):
    """
    打印键长信息
    
    Parameters:
    -----------
    bonds : list
        键长信息列表
    sort_by : str
        排序方式: 'distance', 'element', 'index'
    """
    
    if sort_by == 'distance':
        bonds.sort(key=lambda x: x[4])
    elif sort_by == 'element':
        bonds.sort(key=lambda x: (x[2], x[3], x[4]))
    elif sort_by == 'index':
        bonds.sort(key=lambda x: (x[0], x[1]))
    
    print(f"\n{'原子1':<8} {'原子2':<8} {'元素1':<6} {'元素2':<6} {'键长(Å)':<12}")
    print("-" * 50)
    
    for i, j, elem1, elem2, distance in bonds:
        print(f"{i+1:<8} {j+1:<8} {elem1:<6} {elem2:<6} {distance:<12.4f}")


def save_bond_lengths(bonds, filename, sort_by='distance'):
    """
    保存键长信息到文件
    
    Parameters:
    -----------
    bonds : list
        键长信息列表
    filename : str
        输出文件名
    sort_by : str
        排序方式
    """
    
    if sort_by == 'distance':
        bonds.sort(key=lambda x: x[4])
    elif sort_by == 'element':
        bonds.sort(key=lambda x: (x[2], x[3], x[4]))
    elif sort_by == 'index':
        bonds.sort(key=lambda x: (x[0], x[1]))
    
    with open(filename, 'w') as f:
        f.write("# 键长分析结果\n")
        f.write(f"# 总共找到 {len(bonds)} 个键\n\n")
        f.write(f"{'原子1':<8} {'原子2':<8} {'元素1':<6} {'元素2':<6} {'键长(Å)':<12}\n")
        f.write("-" * 50 + "\n")
        
        for i, j, elem1, elem2, distance in bonds:
            f.write(f"{i+1:<8} {j+1:<8} {elem1:<6} {elem2:<6} {distance:<12.4f}\n")


def analyze_bond_types(bonds):
    """
    分析不同类型的键
    
    Parameters:
    -----------
    bonds : list
        键长信息列表
    
    Returns:
    --------
    bond_types : dict
        不同键类型的统计信息
    """
    
    bond_types = {}
    
    for i, j, elem1, elem2, distance in bonds:
        bond_type = '-'.join(sorted([elem1, elem2]))
        
        if bond_type not in bond_types:
            bond_types[bond_type] = []
        
        bond_types[bond_type].append(distance)
    
    return bond_types


def print_bond_statistics(bond_types):
    """
    打印键长统计信息
    
    Parameters:
    -----------
    bond_types : dict
        不同键类型的距离列表
    """
    
    print(f"\n{'键类型':<12} {'数量':<8} {'最短(Å)':<12} {'最长(Å)':<12} {'平均(Å)':<12} {'标准差(Å)':<12}")
    print("-" * 80)
    
    for bond_type, distances in sorted(bond_types.items()):
        count = len(distances)
        min_dist = min(distances)
        max_dist = max(distances)
        avg_dist = np.mean(distances)
        std_dist = np.std(distances)
        
        print(f"{bond_type:<12} {count:<8} {min_dist:<12.4f} {max_dist:<12.4f} {avg_dist:<12.4f} {std_dist:<12.4f}")


def main():
    parser = argparse.ArgumentParser(
        description="从CIF文件提取键长信息",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
使用示例:
  python extract_bond_lengths.py structure.cif
  
  python extract_bond_lengths.py structure.cif --max-distance 5.0 --min-distance 1.0
  
  python extract_bond_lengths.py structure.cif --output bonds.txt
  
  python extract_bond_lengths.py structure.cif --sort element
  
  python extract_bond_lengths.py structure.cif --stats-only
        """
    )
    
    parser.add_argument("cif_file", help="输入的CIF文件")
    parser.add_argument("--max-distance", type=float, default=4.0,
                       help="最大键长阈值 (Å) (默认: 4.0)")
    parser.add_argument("--min-distance", type=float, default=0.5,
                       help="最小键长阈值 (Å) (默认: 0.5)")
    parser.add_argument("--output", "-o", help="输出文件名 (可选)")
    parser.add_argument("--sort", choices=['distance', 'element', 'index'], 
                       default='distance', help="排序方式 (默认: distance)")
    parser.add_argument("--stats-only", action="store_true",
                       help="只显示统计信息,不显示详细键长")
    
    args = parser.parse_args()
    
    try:
        atoms = read(args.cif_file)
        print(f"成功读取文件: {args.cif_file}")
        print(f"化学式: {atoms.get_chemical_formula()}")
        print(f"原子数: {len(atoms)}")
        print(f"晶胞参数: {atoms.get_cell_lengths_and_angles()}")
    except Exception as e:
        print(f"错误: 无法读取文件 {args.cif_file}")
        print(f"错误信息: {e}")
        return
    
    print(f"\n正在计算键长...")
    print(f"键长范围: {args.min_distance:.2f} - {args.max_distance:.2f} Å")
    
    bonds = get_bond_lengths(atoms, 
                           max_distance=args.max_distance,
                           min_distance=args.min_distance)
    
    print(f"找到 {len(bonds)} 个键")
    
    bond_types = analyze_bond_types(bonds)
    
    print_bond_statistics(bond_types)
    
    if not args.stats_only:
        print_bond_lengths(bonds, sort_by=args.sort)
    
    if args.output:
        save_bond_lengths(bonds, args.output, sort_by=args.sort)
        print(f"\n键长信息已保存到: {args.output}")
        
        stats_file = args.output.replace('.txt', '_stats.txt')
        with open(stats_file, 'w') as f:
            f.write("# 键长统计信息\n\n")
            f.write(f"文件: {args.cif_file}\n")
            f.write(f"化学式: {atoms.get_chemical_formula()}\n")
            f.write(f"原子数: {len(atoms)}\n")
            f.write(f"键长范围: {args.min_distance:.2f} - {args.max_distance:.2f} Å\n")
            f.write(f"总键数: {len(bonds)}\n\n")
            
            f.write(f"{'键类型':<12} {'数量':<8} {'最短(Å)':<12} {'最长(Å)':<12} {'平均(Å)':<12} {'标准差(Å)':<12}\n")
            f.write("-" * 80 + "\n")
            
            for bond_type, distances in sorted(bond_types.items()):
                count = len(distances)
                min_dist = min(distances)
                max_dist = max(distances)
                avg_dist = np.mean(distances)
                std_dist = np.std(distances)
                
                f.write(f"{bond_type:<12} {count:<8} {min_dist:<12.4f} {max_dist:<12.4f} {avg_dist:<12.4f} {std_dist:<12.4f}\n")
        
        print(f"统计信息已保存到: {stats_file}")


if __name__ == "__main__":
    main()
