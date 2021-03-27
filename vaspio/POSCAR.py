__author__ = 'youmu'
__email__ = 'railgun@xjtu.edu.cn'
__date__ = '2021/03/26'
import numpy as np
def skip_n(f, n):
    [f.readline() for i in range(n)]

class POSCAR:
    '''
    Class `POSCAR` is implemented to read, modify and write POSCAR/CONTCAR files.
    To read a POSCAR or CONTCAR, specify its path.
    Make sure to install ase package before modifying and writing POSCAR/CONTCAR files  
    The format of POSCAR/CONTCAR is given at the vaspwiki: https://www.vasp.at/wiki/index.php/POSCAR
    '''
    
    def __init__(self, file=None):

        if file is not None:
            self.read_poscar(file)
            
    def read_poscar(self, file):
        '''
        Take the path of the file which is to be read.
        Note that the extra information such the velocities of ions generated from MD
        simulations will not be read.
        '''
        self.file = file
        try:
            self.read_header()
            self.read_coords()
            self.read_relaxiation_array()
        except FileNotFoundError:
            print('%s file not found, please check the path!'%self.file)
    
    def read_header(self):
        '''
        Read the header of POSCAR file, these following parameters will be loaded:
        comment line
        universal scaling factor
        unit cell
        element types
        number of per element
        Mode of selective calculation
        Mode of coordinates, Cartesian or Direct
        '''
        with open(self.file) as f:
            self.comment = f.readline().strip()  # remove whitespace by strip method.
            self.scaling = float(f.readline().split()[0])
            self.read_cell()
            skip_n(f, 3) # skip over the basis vectors
            self.elements = f.readline().split()
            self.n_per_ion = np.array(f.readline().split(), dtype=int)
            self.n_ions = self.n_per_ion.sum()
            xstr = f.readline()[0]
            if xstr in ['s', 'S']:
                self.selective = True
                self.cart = POSCAR.if_cart(f.readline())
            else:
                self.selective = False
                self.cart = POSCAR.if_cart(xstr)
        self.direct = not self.cart
    def read_cell(self):
        self.cell = np.genfromtxt(self.file, skip_header=2, max_rows=3, dtype=np.float64)
    def read_coords(self):
        '''
        Read atoms' coordinates, return a (n_ions, 3) np.ndarray.
        '''
        self.coords = np.genfromtxt(self.file, skip_header=8+self.selective, max_rows=self.n_ions, usecols=[0, 1, 2], dtype=np.float64)
        
    def read_relaxiation_array(self):
        '''
        Read the states of atoms of being free or fixed, if the selective mode is on.
        Return a boolean array with same shape of ions' coordinates, Ture for free and False for fixed.
        '''
        if self.selective:
            with open(self.file) as f:
                relaxiation_array = []
                skip_n(f, 8+self.selective)
                for i in range(self.n_ions):
                    fix_row = f.readline().split()[3:6]
                    for index, xstr in enumerate(fix_row):
                        if xstr == 'T':
                            fix_row[index] = True
                        elif xstr == 'F':
                            fix_row[index] = False
                        else:
                            raise ValueError("At line %d, F or T need to be set to fix the ion"%(index+8+self.selective))
                    relaxiation_array.append(fix_row)
            self.relaxiation_array = np.array(relaxiation_array, dtype=bool)
            
    def get_index(self, element=None, x=None, y=None, z=None):
        '''
        Get the atom indices, selection conditions can be element type or space range or both.
        An array of index, which starts from 0, of the selected atoms will be returned.
        ***********************************************************************************************
        Input:
        element    :  element type, str
        x          :  squence like, selected range of the first coordinate. For example, x = [0, 0.2]
        y          :  Same as x but the second coordinate
        z          :  Same as x but the third coordinate
        Output:
        An array of the selected atom`s indices, which start from 0.
        ***********************************************************************************************
        
        To generate a mask, try the following code:
        ************************************************
        import numpy as np
        a = POSCAR('./POSCAR')
        indices = a.get_index('Fe')
        mask = np.zeores(a.n_ions, dtype=bool)
        mask[indices] = 1
        ************************************************
        '''
        if element is not None:
            atom_indices = self.get_index_element(element)
            mask = np.zeros(self.n_ions, dtype=bool)
            mask[atom_indices] = 1
        else:
            mask = np.ones(self.n_ions, dtype=bool)
        return np.nonzero(np.logical_and(mask, self.get_index_pos(x, y, z)))[0]
    
    def get_index_element(self, element):
        '''
        Base method of get_index, see more at the description about get_index method.
        Return the atom indices of the selected element, starting from 0. 
        '''
        # element selection is calculated from self.elements and self.n_per_ions
        assert element in self.elements
        indices = [index for index, _element in enumerate(self.elements) if _element == element]
        new_n_per_ions = np.append(np.array([0]), self.n_per_ion)
        n_cumsum = new_n_per_ions.cumsum()
        atom_indices_lst = [np.arange(n_cumsum[index], n_cumsum[index+1]) for index in indices] 
        return np.concatenate(atom_indices_lst) # concatenate method, join a sequence of arrays along an existing axis

        
    
    def get_index_pos(self, x=None, y=None, z=None):
        '''
        Base method of get_index, see more at the description about get_index method.
        Return the atom indices by selecting coordinates range, starting from 0. 
        '''
        if x or y or z is not None:
            return POSCAR.pos_filter(self.coords, x, y, z)  # return a mask
        else: return np.ones(self.n_ions, dtype=bool)
    
    @staticmethod
    def pos_filter(coords, x, y, z):
        '''
        Method to filter atoms by the input coordinates range.
        '''
        pos_range = [x, y, z]
        control = np.array([x, y, z], dtype=bool)  # generate the control array, None will be set to False
        mask = np.zeros_like(coords, dtype=bool)   
        for i in range(control.size):              # True to make logical operations and False to fill the 
            if control[i] == False:                # corresponding mask column to True.
                mask[:,i].fill(1)
            elif control[i] == True:
                _min, _max = np.min(pos_range[i]), np.max(pos_range[i])
                mask[:,i] = np.logical_and(coords[:, i]>=_min, coords[:, i]<=_max)
        return mask.sum(axis=1) > 2   # Sum over the axis 1, value 3 means true. All (x, y, z) must be in the range
    @staticmethod
    def if_cart(astr):
        try:
            float(astr[0])
            raise TypeError('Please input coordination mode!')
        except(ValueError):
            if astr[0] in ['C', 'c','K','k']:
                return True
            else: return False
################# TODO #########################
# instance method, modify and write POSCAR