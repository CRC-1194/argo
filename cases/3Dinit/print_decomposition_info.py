#!/usr/bin/env python
import argparse
import os
import sys

def decomposition_info(vtk_dir, pattern="cell"):
    os.chdir(os.path.join(os.getcwd(), vtk_dir))
    tet_vtks = [tet_vtk for tet_vtk in os.listdir(os.getcwd()) \
                if tet_vtk.startswith(pattern) \
                and tet_vtk.endswith(".vtk")]

    n_interface_cells = len(tet_vtks)

    if n_interface_cells == 0:
        print("Error: no interface cells found in %s. Have the tet-decompositions actually been written to disk?" % vtk_dir)
        sys.exit(1)

    tets_per_cell = []
    accumulated_tets = 0
    for tet_vtk in tet_vtks:
        vtk = open(tet_vtk)
        for line in vtk:
            if "CELL_DATA" in line:
                n_tets = int(line.split()[1])
                tets_per_cell.append(n_tets)
                accumulated_tets = accumulated_tets + n_tets
    
    print("Interface cells:", n_interface_cells)
    print("Average tets per interface cell:", accumulated_tets/float(n_interface_cells))


parser = argparse.ArgumentParser(description='Compute the average number of tets per interface cell for SMCA.')

parser.add_argument('directory', type=str,
                    help='Directory where the tetrahedral cell decompositions are located\n, usually case_dir/VTK')

args = parser.parse_args()

if __name__ == '__main__':
    decomposition_info(args.directory)
