import glob
import os
import sys





extensions = [".cpp", ".hpp"]

PATH = "/home/youngmit/git/mocc/src"

def strip( filename ):
    with open(filename, 'r') as f:
        new = []
        n_white = 0
        for line in f:
            strip_line = line.rstrip()
            # Minus 1 for the newline
            n_white += len(line) - len(strip_line) - 1
            new.append(strip_line)
    with open(filename, 'w') as f:
        [f.write('%s\n' % line) for line in new]

    print "Stripped " + str(n_white) + " characters of whitespace."
    return



for path, dirs, files in os.walk(PATH):
    for f in files:
        file_name, file_extension = os.path.splitext(f)
        if file_extension in extensions:
            while True:
                inp = raw_input("Treat file: " + f + "? [Y/n/q]: ")
                if( inp == "" or inp == "y" or inp == "Y" ):
                    strip(os.path.join(path, f))
                    break
                if( inp == "n" or inp == "N" ):
                    # next
                    continue
                if( inp == "q" ):
                    exit()





