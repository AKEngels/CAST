import os

ENDINGS = (".cc", ".h")

# function to count number of lines of a file
def count_lines(filename):
    with open(filename) as inp:
        lines = inp.readlines()
    return len(lines)

number_of_lines = 0
for root, dirs, files in os.walk("."):       # also search all subdirectories
    for name in files:
        if name.endswith(ENDINGS):           # for files with endings defined above (source files)
            filepath = root + os.sep + name
            number_of_lines += count_lines(filepath) # add the number of lines

print number_of_lines
