### script that searches the lowest path through a PES or a Free Energy Surface
### the surface is given as a csv file (with labels, see 'test.csv' as an example)
### results are a latex document 'pes.tex' which contains a table
### where the lowest energy path is marked blue
### a png file where the energy of the path is plotted
### and a csv file 'path.csv' which contains a counter, 
### the x- and y-coordinate as well as the energy of every point on the path

# USER INPUT
START = "up_left"            # in which corner to start?
FILENAME = "freeEnergy.csv"  # name of the csv file

################################################################################
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# tests if a string can be converted to a float
def is_float(string):
    try:
        float(string)
    except ValueError:
        return False
    else:
        return True

def add_one(var):
    return var+1

def sub_one(var):
    return var-1

# first return value is for x, second for y
def get_move_functions(start):
    if start == "down_left":  
        return add_one, sub_one
    elif start == "up_right": 
        return sub_one, add_one
    elif start == "up_left": 
        return add_one, add_one
    elif start == "down_right":  
        return sub_one, sub_one
    else:
        print "Unvalid direction!"

# function to get starting coordinates
def get_starting_coordinates(start):
    if start == "down_left":  
        y = len(pes)-2 
        x = 1
    elif start == "up_right": 
        y = 1              
        x = len(pes[0])-2
    elif start == "up_left": 
        y = 1             
        x = 1
    elif start == "down_right":  
        y = len(pes)-2             
        x = len(pes[0])-2
    else:
        print "Unvalid direction!"
    return x,y

# function for condition when walk has finished
def break_condition(x,y,start):
    if start == "down_left":  
        if y > 1 or x < len(pes[0])-2:
            return True
    elif start == "up_right": 
        if x > 1 or y < len(pes)-2:
            return True
    elif start == "up_left": 
        if x < len(pes[0])-2 or y < len(pes)-2:
            return True
    elif start == "down_right":  
        if x > 1 or y > 1:
            return True
    else:
        print "Unvalid direction!"
    return False

# read csv file
with open(FILENAME) as pesfile:
    lines = pesfile.readlines()

# save content of csv file in matrix 'pes'
pes = []
for l in lines:
    linelist = l.split(",")
    for i,c in enumerate(linelist):
        if c.find("\n") != -1:
            c = c.replace("\n","")
            linelist[i] = c
    pes.append(linelist)

# get labels
x_labels = pes[0][1:len(pes[0])-1]  # first line
y_labels = []                       # first column
for i, line in enumerate(pes):
    if i != 0 and i != len(pes)-1:
        y_labels.append(line[0])

# find and write path
path_xy = []   # vector to save the position of the elements that make up the path
energies = []  # vector to save energies of path
with open("path.csv","w") as path:
    path.write("step, x, y, value\n")  # write headline
    counter = 1
    x,y = get_starting_coordinates(START)
    move_x, move_y = get_move_functions(START)
    if (pes[y][x]) == "":
        print "Unvalid starting point. No value here."
    
    # move through the PES at the lowest energy path
    while break_condition(x,y,START): 
        value = float(pes[y][x])  # get energy value of current point...
        path.write("{}, {}, {}, {}\n".format(counter, x_labels[x-1], y_labels[y-1], value)) 
        path_xy.append([x,y])  # ...and save position...
        energies.append(value) # ...and energy
        counter += 1
    
        # possible energy value for next step
        value_y = pes[move_y(y)][x]  
        value_x = pes[y][move_x(x)]   

        # if one of those values if empty the other one is the right one
        if value_x == "":
            y = move_y(y)
        elif value_y == "":
            x = move_x(x)

        # if both are valid energies move to the lower one  
        else:
            value_x = float(value_x)
            value_y = float(value_y)

            if value_y < value_x:
                if move_y == sub_one and y == 1:
                    if x != 1 or move_x == add_one:
                        x = move_x(x)
                else:
                    y = move_y(y)
            elif value_x < value_y:
                if move_x == sub_one and x == 1:
                    if y != 1 or move_y == add_one:
                        y = move_y(y)
                else:
                    x = move_x(x)
                    
            elif value_x == value_y:  # if both are identical...
                path.write("{}, , , {}\n".format(counter, value_x))  # ...write any of them to file...
                path_xy.append([x,move_y(y)])  # ...save both coordinates...
                path_xy.append([move_x(x),y])
                energies.append(value_x) # ...save the energy...
                counter += 1
                y = move_y(y)  # ...and go 2 steps further
                x = move_x(x)

    # save last point
    value = float(pes[y][x]) 
    path.write("{}, {}, {}, {}\n".format(counter, x_labels[x-1], y_labels[y-1], value)) 
    path_xy.append([x,y])  
    energies.append(value)
    counter += 1


# plot path
plt.plot(energies)
plt.xlabel("Step")
plt.ylabel("Energy [kcal/mol]")
plt.savefig("path.png")
plt.close()

# write PES as table in LaTeX file
with open("pes.tex","w") as texfile:
    texfile.write("\documentclass[landscape]{article}\n")        # documentclass
    texfile.write("\usepackage[dvipsnames,table]{xcolor}\n")     # package for color
    texfile.write("\usepackage{geometry}\n")                     # package for margins
    texfile.write("\geometry{a0paper, top=25mm, left=0mm, right=0mm, bottom=0mm,")
    texfile.write(" headsep=0mm, footskip=0mm}\n\n") 
    texfile.write(r"\begin{document}"+"\n")  # begin document
    texfile.write(r"\tiny")                  # small text size
    texfile.write(r"\centering")             # table in the center of page
    
    texfile.write(r"\begin{tabular}{c|")  # headline of table...
    for i in range(1,len(pes[0])):
        texfile.write("c")
    texfile.write("}")
    headline = ""
    for column in pes[0]:  # ...from first line of PES
        headline += str(column) + " & "
    texfile.write(headline[:-2] + r"\\" + "\n")
    texfile.write(r"\hline"+"\n")

    for y,line in enumerate(pes):  # write rest of PES in LaTeX file
        if y > 0:
            linestring = ""
            for x,column in enumerate(line):
                if [x,y] in path_xy:  # if x and y coordinate of current point is saved...
                    rounded = str(round(float(column),1))
                    linestring += r"\cellcolor{blue!25}"+rounded+" & "
                else:
                    if is_float(column) and x > 0:
                        linestring += str(round(float(column),1)) + " & "
                    else:
                        linestring += str(column) + " & "
            texfile.write(linestring[:-2] + r"\\" + "\n")

    texfile.write(r"\end{tabular}"+"\n")  # end of table
    texfile.write(r"\end{document}")      # end of document
