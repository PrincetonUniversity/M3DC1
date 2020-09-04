import sys

f1 = open(str(sys.argv[1]), "r")
f2 = open(str(sys.argv[2]), "r")
exitval = 0
tol = 1e-3

cols = ["ntime", "time", "ekin", "gamma_gr", \
        "ekinp", "ekint", "ekin3", \
        "emagp", "emagt", "emag3", \
        "etot"]

while True:
    v1 = f1.readline().split()
    v2 = f2.readline().split()

    if len(v1) is 0:
        if len(v2) is not 0:
            print("Error: files are different lengths")
            exitval = 1
            break
        break
    
    if len(v1) is not len(v2):
        print("Error: lines are different lengths")
        exitval = 1
        break

    for c, x, y in zip(cols, v1, v2):
        if c == "etot":
            continue
        
        if float(x) == 0.0:
            if float(y) != 0.0:
                print("Files differ at time", v1[0])
                exitval = 1
                break
        elif abs((float(x) - float(y)) / float(x)) > tol:
            print("Files differ at time", v1[0])
            print(c, " (base) = ", float(x))
            print(c, " (new) = ", float(y))
            print("Fractional difference = ", \
                  abs((float(x) - float(y)) / float(x)))
            print("Tolerance = ", tol)
            exitval = 1
            break
        
    if exitval is 1:
        break

sys.exit(exitval)
