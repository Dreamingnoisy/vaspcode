import math,matplotlib.pyplot as plt

def distance(a,b,cell): #calculate the distance two atoms under periodic boundary condition
    #Assume the cell is orthogonal
    dx = abs(a[0] - b[0])
    dx = min(dx, abs(cell[0] - dx))

    dy = abs(a[1] - b[1])
    dy = min(dy, abs(cell[1] - dy))

    dz = abs(a[2] - b[2])
    dz = min(dz, abs(cell[2] - dz))

    return math.sqrt(dx * dx + dy * dy + dz * dz)

def getXYZ(str):
    '''
    return the x, y, z coordinates of one line in the form of a tuple
    Example: getXYZ('Fe  0   0.5   0.5') -----> (0,0.5,0.5)
    '''
    XYZ = []
    for value in str.split()[1:]:
        XYZ.append(float(value)) 
    return tuple(XYZ)

def getVolume(z,r,z_top,z_bot):
    volume = 4/3*math.pi*r**3
    if z+r > z_top:
        h = z+r - z_top        
        volume -= math.pi*h**2*(r - h/3)
    if z-r < z_bot:
        h = z_bot - (z - r)
        volume -= math.pi*h**2*(r - h/3)     
    return volume

class Trajectory:
    def __init__(self,filename,cell,z_top,z_bot,*elements,resolution=200,skip=0):
        '''
        filename  : path to the trajectory file, .xyz form required\n
        cell      : lattice length of 3 axis of the cell\n
        z_top and z_bot are set to handle the occasion that a volume between z_top and z_bot \n
        need to be selected.\n
        z_top     : average vertical coordinate for upper interface\n
        z_bot     : average vertical coordinate for lower interface\n
        elements  : select two elements to calculate the g_of_r\n
        resolution: number of points in the final radial distribution function\n
        skip      : determine the step size of reading the trajectory file. Use skip = 0 to read all steps of MD \n
        '''

        self.filename = filename
        self.z_top,self.z_bot = z_top,z_bot
        self.elements = elements
        self.resolution = resolution
        self.skip = skip
        self.cell =cell   # A, B, C of the lattice
        self.read_str()
        self.rho_0_cal()

    def __str__(self):
        str0 = 'A Trajectory class object reading from %s\n'%self.filename
        str1 = 'The atom list is %s\n'%str(self.atom_list)
        str2 = '%s element(s) has been selected\n'%str(self.elements)
        str3 = '%s images have been read.'%self.n_image
        return str0 + str1 + str2+str3

    def read_str(self):
        self.select_element()

        f = open(self.filename, 'r')
        self.data = f.readlines()
        f.close()

        self.n_atoms = int(self.data[0].split()[0])  #read the total number of atoms
        l_image = self.n_atoms + 2            #length of an image
        n_steps = len(self.data) // l_image   # read the number of steps
        self.get_atom_list()                  #generate the list of atoms

        self.coordinates = [] 
        for step in range(0,n_steps,self.skip+1):
            coords_image = self.readimage(2+step*l_image)
            self.coordinates.append(tuple(zip(self.atom_list,coords_image)))  
        self.n_image = len(self.coordinates)  #get the number of images read from self.data
            
    def get_atom_list(self):
        '''
        Extract the list of atoms from the data
        returns a list
        '''
        self.atom_list = []
        for line in self.data[2:self.n_atoms+2]:
            self.atom_list.append(line.split()[0])

    def readimage(self,beg):
        '''
        read the XYZ coordinates of one image from the 'beg' line
        return a list of coordinates
        '''
        coords_image = []
        for line in self.data[beg:beg+self.n_atoms]:  
            coords_image.append(getXYZ(line))
        return coords_image
    
    def select_element(self):
        '''
        Purse the self.elements, set some default values and raise a warning.
        '''

        if len(self.elements) > 2:
            print("Too many elements are selected! Two maximum elements are required.")
            raise UserWarning(self.elements)
        if len(self.elements) == 2:
            print('Two element types are selected as the element pair: %s , %s'%(self.elements[0],self.elements[1]))

        if len(self.elements) == 0:
            print('The element type has not been specified.\nChoose %s %s as the default element pair'%(self.atom_list[0],self.atom_list[0]))
            self.elements = (self.atom_list[0],self.atom_list[0])
        
        if len(self.elements) == 1:
            print('Only one element type is selected %s:'%self.elements[0])
            print('Set %s %s as the element pair:'%(self.elements[0],self.elements[0]))
            self.elements = (self.elements[0],self.elements[0])
        
    def rho_0_cal(self):
        '''
        Calculate the average density of each selected atom

        rho_0 = Volume / num_of_selected_pairs
        for same elements : Npair = N1*(N1-1)
        for different element: Npair = N1*N2
        '''
        n_elements0,n_elements1 = 0,0
        for image in self.coordinates:
            for atom in image:
                if atom[0] == self.elements[0]:
                    if  self.z_bot  <  atom[1][2] < self.z_top:
                        n_elements0 += 1
                if atom[0] == self.elements[1]:
                    if  self.z_bot  <  atom[1][2] < self.z_top:
                        n_elements1 += 1
        volume = self.cell[0]*self.cell[1]*(self.z_top - self.z_bot)*self.n_image

        if self.elements[0] == self.elements[1]:    
            self.rho_0 = (n_elements0*(n_elements1-1*self.n_image))/volume
        else:
            self.rho_0 = (n_elements0*n_elements1)/volume   
        # the rho_0 are amplified by a factor of self.n_image because of the multiplication between N1*(N1-1) or N1 * N2
        self.rho_0 = self.rho_0/self.n_image
    def getelements_data(self,image):
        # extract the coordinates of the selected elements from each image
        coords_atom0,coords_atom1 =[],[]
        for atom in image:
            if atom[0] == self.elements[0]:
                if  self.z_bot  <  atom[1][2] < self.z_top:
                    coords_atom0.append(atom)
            if atom[0] == self.elements[1]:
                if  self.z_bot  <  atom[1][2] < self.z_top:
                    coords_atom1.append(atom)
        return coords_atom0,coords_atom1

    def rdf_cal(self):
        self.r_max = min(self.cell[0],self.cell[1]) / 2.0
        self.dr = self.r_max / self.resolution
        g_of_r_perimage = []

        #loop over all the images in the self.coordinates
        #for each image, the values g_of_r at each point are stored in g_of_r_perimage.
        #the g_of_r of each image can be exported if necessary
        for image in self.coordinates:
            coords_atom0,coords_atom1 = self.getelements_data(image)
            dV_per_atom0,g_of_r = [],[0]*self.resolution
            for atom0 in coords_atom0:
                dV_per_r = []
                for i in range(self.resolution):
                    r = i*self.dr
                    dV = getVolume(atom0[1][2],r+self.dr,self.z_top,self.z_bot) - getVolume(atom0[1][2],r,self.z_top,self.z_bot)
                    dV_per_r.append(dV)
                dV_per_atom0.append(dV_per_r)

                for atom1 in coords_atom1:
                    if atom0 != atom1:
                        dist = distance(atom0[1],atom1[1],self.cell)
                        index = int(dist/self.dr)
                        if 0< index < self.resolution:   
                            g_of_r[index] +=1/dV_per_r[index]/self.rho_0
            g_of_r_perimage.append(g_of_r)

        self.g_of_r =[0]*self.resolution
        for image in g_of_r_perimage:    #ensomble average, sum the g_of_r_perimage and normalize by the factor of
            for i in range(len(image)):  #the length of the images read from self.coordinates 
                self.g_of_r[i] += image[i]/self.n_image

    def rdf_plt(self,outfile=None):
        Rr = []
        for i in range(self.resolution):
            Rr.append(i*self.dr)

        plt.xlabel('r (Å)')
        plt.ylabel('g$_{%s-%s}$(r)'%(self.elements[0],self.elements[1]))    
        plt.plot(Rr,self.g_of_r)
        if outfile is not None:
            plt.savefig(outfile,dpi=300,bbox_inches='tight')
        plt.show()

if __name__ == '__main__':
    cell = 10,10,21
    FeCrNi = Trajectory('movie.xyz',cell,21.0,0.0,'Fe','Fe')
    FeCrNi.rdf_cal()
    FeCrNi.rdf_plt('rdf.png')
    
