STRUCTURE = "19_21.arc"
MD_FILE = "output_MD_SNAP.arc"

FIX_ATOMS = "1-78,237-461,594-912,994-1004,1042-1143,1233-1324,1396-1437,1462-1679,1690-1727,1807,1809-1815,1819-1824,1828-1845,1849-1851,1855-1866,1870-1881,1888-1896,1900-1923,1927-1932,1936-2010,2014-2034,2038-2043,2047-2049,2056-2058,2062-2064,2068-2070,2077-2085,2089-2103,2107-2139,2143-2151,2155-2169,2173-2184,2188-2205,2209-2226,2233-2250,2257-2259,2266-2280,2284-2292,2302-2337,2341-2379,2383-2388,2392-2394,2398-2403,2407-2418,2422-2424,2428-2436,2440-2493,2497-2520,2527-2538,2545-2547,2554-2583,2587-2598,2602-2613,2620-2628,2632-2637,2641-2670,2674-2679,2683-2688,2698-2727,2734-2748,2752-2769,2773-2817,2821-2841,2854-2856,2860-2892,2896-2901,2905-2916,2920-2925,2929-2940,2944-2949,2959-2979,2983-2988,2992-3006,3010-3015,3019-3042,3052-3060,3064-3066,3073-3075,3079-3111,3115-3123,3127-3129,3139-3153,3157-3159,3163-3165,3169-3174,3181-3198,3205-3207,3211-3213,3217-3222,3226-3246,3250-3267,3271-3309,3316-3339,3343-3348,3352-3354,3358-3369"
NUMBER_OF_ATOMS = 3369

import math

def parse_range(string):
    rangelist = string.split(",")
    ranges = []
    for r in rangelist:
        rlist = r.split("-")
        if len(rlist) == 1:
            ranges += [int(rlist[0])]
        else:
            start = int(rlist[0])
            stop = int(rlist[1])
            ranges += range(start, stop+1)
    return ranges

def get_dist(line1, line2):
    x1 = float(line1.split()[2])
    y1 = float(line1.split()[3])
    z1 = float(line1.split()[4])
    x2 = float(line2.split()[2])
    y2 = float(line2.split()[3])
    z2 = float(line2.split()[4])
    res = math.sqrt((x1-x2)*(x1-x2) + (y1-y2)*(y1-y2) + (z1-z2)*(z1-z2))
    return res

def calc_rmsd(structure1, structure2, ACT_ATOMS):

    if len(structure1) != len(structure2):
        print "Those structures can't be compared"
        return 0

    else:
        delta_square = []
        for i, line in enumerate(structure1):
            if i in ACT_ATOMS:
                delta_square.append(get_dist(line, structure2[i])*get_dist(line, structure2[i]))
        rmsd = math.sqrt(sum(delta_square)/len(ACT_ATOMS))
        return rmsd

def get_act_atoms():
    fixed_atoms = parse_range(FIX_ATOMS)
    act_atoms = []
    for a in TOTAL_ATOMS:
        if a not in fixed_atoms:
            act_atoms.append(a)
    return act_atoms

#TOTAL_ATOMS = range(1,NUMBER_OF_ATOMS+1)  # numbering starting with 1
ACT_ATOMS = parse_range("216, 217, 218, 219, 220, 1009, 1010, 1011, 1012, 1013, 1014, 1015, 1016, 1017, 1018, 1019, 1020, 1021, 1022, 1208, 1209, 1210, 1211, 1212, 1213, 1214, 1215, 1216, 1217, 1218, 1219, 1220, 1365, 1366, 1367, 1368, 1369, 152-155,1196-1207,1746-1750,1758-1768,1779-1786,1804-1806,86-93,146-151,156-161,523-535,1178-1186,1190-1195,1351-1358,1376-1393")

with open(STRUCTURE) as strucfile:
    ORIGINAL_STRUCTURE = strucfile.readlines()

with open(MD_FILE) as mdfile:
    lines = mdfile.readlines()

NUMBER_OF_STRUCTURES = len(lines) / (NUMBER_OF_ATOMS + 1)

rmsds = []
for struc in range(NUMBER_OF_STRUCTURES):
    structure = lines[struc*(NUMBER_OF_ATOMS+1) : (struc+1)*(NUMBER_OF_ATOMS+1)]
    rmsds.append(calc_rmsd(ORIGINAL_STRUCTURE, structure, ACT_ATOMS))

with open("rmsd_qm_small.csv","w") as rmsdfile:
    for r in rmsds:
        rmsdfile.write(str(r)+"\n")
          
            
