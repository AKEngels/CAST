
### plotting distributions and scatterplot
### meanwhile looking if all files are okay

### can only be used after copying all umbrella.txt files into one folder
### i.e. script 'umbrella2d_analysis.py'

import glob
import math
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# USER INPUT
WINDOWS_X = range(-180, 185, 5)       # values for first restraint (X)
WINDOWS_Y = range(-180, 185, 5)       # values for first restraint (Y)

broken_files = []

for VALUE in WINDOWS_X:
    print "searching distributions with X =", VALUE
    xi_lists = []
    for filename in glob.glob("umbrella_{}_*.txt".format(VALUE)):
        print "looking at file", filename
        with open(filename) as umbrella_file:
            lines = umbrella_file.readlines()
        xi = []
        for i,line in enumerate(lines):
            if float(line.split()[1]) == 0 and float(line.split()[2]) == 0 and filename not in broken_files:
                broken_files.append(filename) # both coordinates are exactly 0 (should not happen accidently)
                break
            xi.append(float(line.split()[2]))
        xi_lists.append(xi)

    plt.rcParams['figure.figsize'] = [11,8]  # 11x8 inches is next to DinA4
    number = int(math.sqrt(len(xi_lists[0])))
    n, bins, patches = plt.hist(xi_lists, number, histtype='step')
    plt.savefig("distribution_X{}.png".format(VALUE), dpi=100)
    plt.close()
    
for VALUE in WINDOWS_Y:
    print "searching distributions with Y =", VALUE
    xi_lists = []
    for filename in glob.glob("umbrella_*_{}.txt".format(VALUE)):
        print "looking at file", filename
        with open(filename) as umbrella_file:
            lines = umbrella_file.readlines()
        xi = []
        for i,line in enumerate(lines):
            if float(line.split()[1]) == 0 and float(line.split()[2]) == 0 and filename not in broken_files:
                broken_files.append(filename) # both coordinates are exactly 0 (should not happen accidently)
                break
            xi.append(float(line.split()[1]))
        xi_lists.append(xi)

    plt.rcParams['figure.figsize'] = [11,8]  # 11x8 inches is next to DinA4

    number = int(math.sqrt(len(xi_lists[0])*len(xi_lists)))

    n, bins, patches = plt.hist(xi_lists, number, histtype='step')
    plt.savefig("distribution_Y{}.png".format(VALUE), dpi=100)
    plt.close()


points = []
print "starting scatterplots"
for filename in glob.glob("umbrella_*.txt"):
    print "looking at file", filename
    with open(filename) as umbrella_file:
        lines = umbrella_file.readlines()
        
    points_x = []
    points_y = []
    for line in lines:
        linelist = line.split()
        points_x.append(float(linelist[1]))
        points_y.append(float(linelist[2]))

    plt.scatter(points_x,points_y)
    plt.savefig("{}.png".format(filename))
    plt.close()
    
    points.append([points_x, points_y, filename])

for i,p in enumerate(points):
    plt.scatter(p[0], p[1], s=0.2, label=p[2])
plt.legend(loc='upper right', fontsize=5)
plt.savefig('scatter.png'.format(filename), dpi=100)


with open("broken_files.txt","w") as outfile:
    outfile.write(str(broken_files))

