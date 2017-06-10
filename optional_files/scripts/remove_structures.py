"""removes a number of structures from the beginning of an MD snapshot file
useful if the file is too big to open it in VMD or if you want to upload it somewhere"""

NAME_OF_THE_SNAPSHOT_FILE = "water_crash.arc"
NUMBER_OF_REMOVED_STRUCTURES = 900
NUMBER_OF_ATOMS = 1261

with open(NAME_OF_THE_SNAPSHOT_FILE, "r") as output:
    counter = 0
    lines = []
    while True:  # loop over all lines
        try:
            x = output.readline() # read line
            if counter >= NUMBER_OF_REMOVED_STRUCTURES * (NUMBER_OF_ATOMS + 1):
                lines.append(x)  # save line after the structures that should be removed
        except:
            break  # stop looping when there are no more lines
        counter = counter + 1
        
       
with open("MD_snap_new.arc", "w") as new: # write saved lines in new file
    for i in lines:
        new.write(i)
    

    
                                                  
                                                  
            

    
