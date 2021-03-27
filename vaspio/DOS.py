__author__ = 'youmu'
__email__ = 'railgun@xjtu.edu.cn'
__date__ = '2021/03/26'

import numpy as np
class DOS():
    def __init__(self, file):
        self.file = file
        self.load_DOSCAR()
        self.read_Header()
        self.get_energies()
        
    def load_DOSCAR(self):
        self._dos = open(self.file, 'rt') 
        
    def read_Header(self):
        self._dos.seek(0)  # shift to 0 byte
        self.n_ions, self.pdos =  np.array(self._dos.readline().split(), dtype=int)[[1, 2]]
        
        [self._dos.readline() for i in range(4)] # skip to line 6
        self.e_max, self.e_min, nedos, self.e_fermi = np.array(self._dos.readline().split()[:-1], dtype=float)
        self.nedos = np.int(nedos)
        
        self.spin = 1 if len(self._dos.readline().split())==5 else 0
        self._dos.close()
        
    def get_energies(self):
        # get the energies and shift it to the Fermi energy
        self.energies = np.genfromtxt(self.file, skip_header=6, max_rows=self.nedos, usecols=0, ) - self.e_fermi
        
    def get_tdos(self):
        if self.spin:
            self.tdos = np.genfromtxt(self.file, skip_header=6, max_rows=self.nedos, usecols=(1, 2)).T
            self.tdos = self.tdos*np.array([1,-1])[:,None]
        if not self.spin:
            self.tods = np.genfromtxt(self.file, skip_header=6, max_rows=self.nedos, usecols=1).T
            
    def get_atom_dos(self, index):
        '''Return the dos of an atom with the given index, which starts from 0'''
        assert self.pdos, "Non partial dos is not available"
        if self.spin:
            a = np.genfromtxt(self.file, skip_header=5+(self.nedos + 1)*(index+1)+1, max_rows=self.nedos)[:,1:]
            return np.array([a[:,::2], -a[:,1::2]])
        else:
            return np.genfromtxt(self.file, skip_header=5+(self.nedos + 1)*(index+1)+1, max_rows=self.nedos)[:,1:]
    
    def sum_atoms_dos(self, indices):
        '''Sum the DOS over all speficifed atoms'''
        dos_sum = np.zeros_like(self.get_atom_dos(indices[0]))
        for index in indices:
            dos_sum= dos_sum + self.get_atom_dos(index)
        return dos_sum
    
    def get_all_pdos(self):
        '''Sum over all atoms' pdos'''
        return self.sum_atoms_dos(np.arange(self.n_ions))
# define orbital dict
orbital_dic = {}
orbital_str = ['s', 'py','pz','px','dxy', 'dyz', 'dz2-r2', 'dxz', 'dx2-y2','fy3x2','fxyz','fyz2', 'fz3', 'fxz2', 'fzx2', 'fx3']
[orbital_dic.setdefault(orbital_str[i], i)  for i in range(len(orbital_str))]

def get_orbitals(orbitals):
    '''
    Return an orbital list corresponding to the input orbitals
    *********************************************************
    Input:
    orbitals: squence like, which contains the aimed orbitals
    *********************************************************
    Example:
    orb_lst = get_orbitals(['s','px'])
    Return:
    [0, 1]
    '''
    orbital_lst = []
    
    for key in orbitals:
        if key not in orbital_dic.keys():
            raise ValueError("%s is not a orbital"%key)
    
    [orbital_lst.append(orbital_dic.get(key)) for key in orbitals]
    return orbital_lst

def get_orbitals_dos(dos, orbitals):
    '''
    Return an ndarray corresponding to the input orbitals
    *********************************************************
    Input:
    dos: ndarray generated from the class methods of class DOS
    orbitals: squence like, which contains the aimed orbitals
    *********************************************************
    Example:
    orb_lst = get_orbital_dos(c,['s','px'])
    Return:
    An ndarray whose first and second columns are the dos of s and 
    p orbitals, respectively.
    '''
    orbital_lst = get_orbitals(orbitals)
    if dos.ndim == 2:
        return dos[:,orbital_lst[:dos.shape[-1]]]
    if dos.ndim == 3:
        return dos[:,:,orbital_lst[:dos.shape[-1]]]
    else:
        raise ValueError('Incorrect dimension detected, the ndim of the dos array should be 2 or 3.\nNow it is %d.'%dos.ndim)
        
def sum_orbitals_dos(dos,orbitals):
    '''
    Sum dos array over all given orbitals
    '''
    if type(orbitals) == str:
        if orbitals == 'S':
            orbitals = ['s']
        elif orbitals == 'P':
            orbitals = ['py','pz','px']
        elif orbitals == 'D':
            orbitals = ['dxy', 'dyz', 'dz2-r2', 'dxz', 'dx2-y2']
        elif orbitals == 'F':
            orbitals = ['fy3x2','fxyz','fyz2', 'fz3', 'fxz2', 'fzx2', 'fx3']
        elif orbitals == 'ALL':
            orbitals = orbital_str
        else:
            raise ValueError('Wrong str received!')
    return get_orbitals_dos(dos, orbitals).sum(axis=-1)       
