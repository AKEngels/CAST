import shutil

def remove_spaces(string):
    """removes unnecessary spaces from a string

    input: string
    return: string without the unnecessary spaces"""

    newstring = ""
    for i, a in enumerate(string):
        if i == 0 and a == " ":  # spaces at the beginning
            pass
        elif a == " " and i == len(string) - 1:  # spaces at the end
            return remove_spaces(newstring)
        elif a == " " and string[i - 1] == " ":  # two spaces back to back
            pass
        else:
            newstring = newstring + a
    return newstring

with open("path.csv") as pathfile:
    lines = pathfile.readlines()

for i,line in enumerate(lines):
    if i != 0:
        linelist = line.split(',')
        x = remove_spaces(linelist[1])
        y = remove_spaces(linelist[2])
        shutil.copy("umbrella_{}_{}.txt".format(x,y), "path/umbrella_{}_{}.txt".format(x,y))
