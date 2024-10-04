from time import gmtime, strftime
import matplotlib.pyplot as plt

def get_first_vx(filename):
    points = []
    with open(filename) as file:
        lines = file.readlines()
        for i in range(3,len(lines)):
            line = lines[i]
            if line == "Vz\n" or line == "Vy\n":
                break
            points.append(float(line))
    return points

def get_nth_vx(filename, n):
    points = []
    with open(filename) as file:
        lines = file.readlines()
        count = 0
        Vx = False
        for i in range(2,len(lines)):
            line = lines[i]
            if line == "Vz\n" or line == "Vy\n" or line == "\n":
                Vx = False
            elif line == "Vx\n":
                Vx = True
                count += 1
            elif Vx and count == n:
                points.append(float(line))
    return points

def get_nth_vy(filename, n):
    points = []
    with open(filename) as file:
        lines = file.readlines()
        count = 0
        Vy = False
        for i in range(2,len(lines)):
            line = lines[i]
            if line == "Vx\n" or line == "Vz\n" or line == "\n":
                Vy = False
            elif line == "Vy\n":
                Vy = True
                count += 1
            elif Vy and count == n:
                points.append(float(line))
    return points

def get_nth_vz(filename, n):
    points = []
    with open(filename) as file:
        lines = file.readlines()
        count = 0
        Vz = False
        for i in range(2,len(lines)):
            line = lines[i]
            if line == "Vx\n" or line == "Vy\n" or line == "\n":
                Vz = False
            elif line == "Vz\n":
                Vz = True
                count += 1
            elif Vz and count == n:
                points.append(float(line))
    return points


title = "Sequential+parallel__expanded_kernel,HALO=4"
#title = "Parallel__expanded_kernel,HALO=4"
#title = "Sequential__expanded_kernel,HALO=4"
plt.title(title)
for i in range(1,9):
    plt.subplot(2,4,i)
    points1 = get_nth_vx("receivers.csv", i)
    points2 = get_nth_vx("receivers_mpi.csv", i)
    sub_title = "Vx="+str(i)
    plt.title(sub_title)
    plt.plot(points1, color= "b")
    plt.plot(points2, color= "r")
#plt.subplot(2,4,1)
#plt.legend(["Sequential", "Parallel, n=2"])

time = strftime("%Y-%m-%d-%H-%M-%S", gmtime())
title = title + " all Vx"
#plt.savefig("figures/visualized/"+title+"__"+time)
plt.show()
            