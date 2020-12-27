import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.interpolate import griddata
def image_point_generator(coord):
    '''
    coord is of direct coordnates
    '''
    x,y,z = coord[0],coord[1],coord[2]
    x1 = x + 1 - (x //0.5) * 2
    y1 = y + 1 - (y //0.5) * 2
    image_coords = []
    image_coords.append((x,y,z))
    image_coords.append((x1,y,z))
    image_coords.append((x,y1,z))
    image_coords.append((x1,y1,z))
    return image_coords

def vector_transform(v,cell):
    v = np.linalg.inv(cell).dot(v)
    v = np.where(np.abs(v) > 0.5, v - np.array([1,1,1])*v/np.abs(v), v)
    v = cell.dot(v)    
    return v
def distance(p1,p2,cell):
    '''
    calculate the distance between two atoms under pbc.
    p1&p2      : points coordinates
    cell       : a numpy array [[a1,b1,c1],
                                [a2,b2,c2],
                                [a3,b3,c3]]
                or [[a1,a2,a3],[b1,b2,b3],[c1,c2,c3]].T
    '''
    delta = p1 - p2
    delta = vector_transform(delta,cell)
    d     = cell.dot(delta)
    return np.sqrt((d**2).sum(axis=-1))

def angle(v1,v2,cell =None ,mode = 'none'):
    '''
    calculate the angle of two vectors under pbc
    v1         : 1*3 numpy array, coords difference between two atoms 
    v2         : 1*3 numpy array
    cell       : same with cell in function distance
    mode       : 'none','single' or 'both'.
                 Need input cell when 'single' or 'both' is turn on.
                 When mode = 'single', only the v1 will be transformed under pbc. 
                 When mode is 'both',both v1 and v2 will be transformed under pbc
                 When mode is 'none',none v1 or v2 will be transformed under pbc
    '''
    if mode == 'both':
        v1 = vector_transform(v1,cell)
        v2 = vector_transform(v2,cell)
    if mode == 'single':
        v1 =vector_transform(v1,cell)        
    if mode == 'none':
        pass          
    cosTheta = v1.dot(v2)/(np.linalg.norm(v1)*np.linalg.norm(v2))
    theta = np.arccos(cosTheta) 
    theta = theta /np.pi * 180
    return theta

#time average operator: <>
#https://chem.libretexts.org/Bookshelves/Physical_and_Theoretical_Chemistry_Textbook_Maps/Book%3A_Time_Dependent_Quantum_Mechanics_and_Spectroscopy_(Tokmakoff)/10%3A_Time-Correlation_Functions/10.02%3A_Correlation_Function_from_a_Discrete_Trajectory
#For detail, see the link above
#Time correlation function of bonds number
def C_of_b_b(b,n):
    
    N = len(b)
    C_of_n = 0
    for i in range(N-n):
        C_of_n += b[i]*b[i+n]
    return C_of_n / (N-n)

def C_of_b_b_2(b,n):
    N = len(b)    
    C_of_n = 0
    for i in range(N//2):
        C_of_n += b[i]*b[n+i]
    return C_of_n / (N//2)

# def C_of_db_db(b,n,b_ave):
    
#     N = len(b)
#     C_of_n = 0
#     for i in range(N-n):
#         C_of_n += (b[i]-b_ave)*(b[i+n]-b_ave)
#     return C_of_n / (N-n)

class Atom:
    def __init__(self,index,element,coord):
        self.index = index
        self.element = element
        self.coord   = coord

class Molecule:
    def __init__(self,name,*atom_index):
        '''
        initialize a molecule by specifing its atom indexes
        '''
        self.name       = name
        self.atom_index = atom_index
        self.atom_dic   = {}
    def atom_add(self,atom):
        self.atom_dic[atom.index] = atom
    def getAtom(self,index):
        return None if index == None else self.atom_dic[index]
class Frame:
    def __init__(self,frame_num,atom_lt):
        self.frame_num = frame_num
        self.atom_lt   = atom_lt
    def getAtom(self,index):
        return None if index == None else self.atom_lt[index-1]
        
    def getAtoms(self,indexes,element):
        atom_lt = []
        [atom_lt.append(self.atom_lt[i-1]) for i in indexes]
        [atom_lt.append(atom) for atom in self.atom_lt if atom.element in element]
        return atom_lt

class Trajectory:
    '''
    Read trajectory from .xyz file and
    do some calculations
    '''
    def __init__(self,input_file,cell):
        self.input_file = input_file
        self.cell       = cell
    def openfile(self):
        f = open(self.input_file,'r')
        #read the total number of atoms
        self.n_atoms = int(f.readline().split()[0])   
        l_image = self.n_atoms + 2             
        f.close()

        self.frame_list = []
        atom_list = []
        for i , l in enumerate(open(self.input_file,'r')):
            ii = i % l_image
            if ii == 0 and i > 0:
                self.frame_list.append(Frame(ii//l_image-1,atom_list))
                atom_list = []
            if ii < 2:
                continue
            entry = l.split()
            atom_list.append(Atom(index = ii -1,element = entry[0],coord=np.array(entry[1:],dtype=float)))
        self.frame_list.append(Frame(ii//l_image-1,atom_list))
        self.frame_num = len(self.frame_list)
    def tcf(self,cutoff,*elements):
        counts_of_frame = []
        for frame in self.frame_list:
            atoms_counter = {}
            pairs_counter = 0
            for atom in frame.atom_lt:
                if atom.element in elements:
                    try:
                        atoms_counter[atom.element].append(atom)
                    except KeyError:
                        atoms_counter[atom.element] = [atom]
            for atom_0 in atoms_counter[elements[0]]:
                for atom_1 in atoms_counter[elements[1]]:
                    d = distance(atom_0.coord,atom_1.coord,self.cell)
                    if d < cutoff:
                        pairs_counter += 1
            counts_of_frame.append(pairs_counter)
            frame.pairs = {elements:pairs_counter}

        #calculating <b> and delat_b(t)    
        b_of_t = np.array(counts_of_frame)
        
        C_values = []   
        for i in range(len(b_of_t)//2):
            C_values.append(C_of_b_b_2(b_of_t,i))
            #C_values.append(C_of_b_b(b_of_t,i))
        C_of_t = np.array(C_values)
        C_of_t = C_of_t / C_of_t[0]
        return C_of_t   

    def tcf_plt(self,ax,t,y):
        # yy = []
        # sumnumber = 50
        # for i ,value in enumerate(y[sumnumber//2:-sumnumber//2]):
        #     y_sum = 0
        #     for j in range(sumnumber):
        #         y_sum += y[i+j] / sumnumber
        #     yy.append(y_sum)    
        ax.plot(t,y)
        
    # !!This function takes a function as a parameter!!
    def orientation(self,molecule,config_f,vnorm): 
        '''
        calculate the orientation of a type of molecules
        molecule: Molecule objects
        config_f:configuration function which takes a molecule instance
                and return a vector
        vnorm  : A vector that used for angle calculation
        '''
        angle_lt = []
        for frame in self.frame_list:
            for index in molecule.atom_index:
                molecule.atom_add(frame.getAtom(index))
            v1 = config_f(molecule,self.cell)
            v2 = vnorm
            an = angle(v1,v2,mode='None')
            angle_lt.append(an)
        return np.array(angle_lt)
    #calculate atomic density maps/contour 
    def adm(self,atom_index):
        coords_lt = []
        for frame in self.frame_list[:]:
            coord = frame.atom_lt[atom_index-1].coord
            coord = np.linalg.inv(self.cell).dot(coord)
            coord = np.where(coord < 0, coord + np.array([1,1,1]), coord)
            image_coords = image_point_generator(coord)
            for image_coord in image_coords:
                coord = self.cell.dot(image_coord) 
                coords_lt.append(coord)
        coords = np.array(coords_lt)
        coords = np.array(sorted(coords,key = lambda coord:coord[0]))
        coords = coords[:]
        points = coords[:,0:2]
        z = coords[:,2]
        xmin ,xmax = np.min(points[:,0]),np.max(points[:,0])
        ymin ,ymax = np.min(points[:,1]),np.max(points[:,1])
        x_grid , y_grid = np.mgrid[xmin:xmax:70j,ymin:ymax:70j]

        z_grid = griddata(points, z, (x_grid, y_grid), method='linear')
        
        return x_grid,y_grid,z_grid,coords

def adm_plt(ax,x_grid,y_grid,z_grid):
    plt.contour(x_grid,y_grid,z_grid)
##These orientation functions should be written for one direction of one molecule 
##in a special system. Each time you want to calculate a new direction, you need a
##configuration function.  
def co2_orien_configuration_1(molecule,cell):
    O66 = molecule.getAtom(66)
    O65 = molecule.getAtom(65)
    v1  = vector_transform(O66.coord - O65.coord,cell)
    return v1

def co2_orien_configuration_2(molecule,cell):
    O65 = molecule.getAtom(65)
    O66 = molecule.getAtom(66)
    C67 = molecule.getAtom(67)
    v1  = 0.5*vector_transform(O66.coord - O65.coord,cell)
    O_midpoint = O65.coord + v1
    v2  = vector_transform(O_midpoint - C67.coord,cell)
    return v2

def co2_orien_configuration_3(molecule,cell):
    O65 = molecule.getAtom(65)
    O66 = molecule.getAtom(66)
    C67 = molecule.getAtom(67)
    v1  = vector_transform(O66.coord - C67.coord,cell)
    return v1
    

def main():
    cell = np.array([[11.3727998699999997 ,0,0],
                    [0,8.0417995,0.000000],
                    [0.000000,0.000000,21.173721]]).T

    t91 = Trajectory(r'T91-773-Cr-Cr-bri\newpos.xyz',cell)
    t91.openfile()
    fig , ax= plt.subplots()
    #y = t91.tcf(3,'C','Fe')
    #t = np.linspace(0,20,200)
    #t91.tcf_plt(t,y)
    #####orentation plot#####
    #co2 = Molecule('co2',65,66,67)
    #vnorm = np.array([0,0,1])
    #angle_lt = t91.orientation(co2,co2_orien_configuration_3,vnorm)
    #####plt angle hist nad hist kde######
    #ax.hist(angle_lt,bins = 180,density = True)
    #sns.kdeplot(angle_lt,shade=True)

    x,y,z,coords = t91.adm(67)
    adm_plt(ax,x,y,z)
    ax.scatter(coords[:,0],coords[:,1],c = coords[:,2])
    plt.colorbar()
    plt.show()
if __name__=="__main__":
    main()
