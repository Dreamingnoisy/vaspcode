from trajectory import Trajectory
if __name__ == '__main__':
    cell = 10,10,21
    FeCrNi = Trajectory('movie.xyz',cell,21.0,0.0,'Fe','Fe')
    print(FeCrNi.n_image)
    FeCrNi.rdf_cal()
    FeCrNi.rdf_plt('rdf.png')