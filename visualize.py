from time import gmtime, strftime
import matplotlib.pyplot as plt

def make_plot(filename):
    points = []
    with open(filename) as file:
        lines = file.readlines()
        for i in range(3,len(lines)):
            line = lines[i]
            if line == "Vz\n" or line == "Vy\n":
                break
            points.append(float(line))
    return points



    
    

points1 = make_plot("receivers_HALO=4.csv")
points2 = make_plot("receivers_mpi.csv")
title = "Sequential+parallel__expanded_kernel,HALO=4"
plt.title(title)
plt.plot(points1, color= "b")
plt.plot(points2, color= "r")
plt.legend(["Sequential", "Parallel"])
time = strftime("%Y-%m-%d-%H-%M-%S", gmtime())
plt.savefig("figures/"+title+"__"+time)
plt.show()
            