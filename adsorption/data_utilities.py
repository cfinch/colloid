#!/usr/bin/env python

# --------------- Functions for working with CSV data ----------------
def read_adsorption_data_CSV(textFileName):
    """Imports adsorption data stored in a .CSV file.  
    Data is in three columns:
    Time(sec),Mean Surface Conc(mol/m^2),Std Dev
    First row contains the column labels
    Subsequent rows contain data
    """
    ## Opening the file
    textFile = open(textFileName, 'r')

    # First line--headers (discard)
    words = textFile.readline().split()

    time = []
    conc = []
    stdev = []

    # Data lines
    isData = True
    while isData:
        words = textFile.readline().strip().split(',')
        if words[0] == '':
            isData = False
        else:
            time.append(float(words[0]))
            conc.append(float(words[1]))
            stdev.append(float(words[2]))
        ####
    ####

    textFile.close()
    return time, conc, stdev
####

def write_adsorption_data_CSV(textFileName, times, mean, std, column_labels):
    """Exports adsorption data to a .CSV file.  
    Data is in three columns:
    Time(sec),Mean Surface Conc(mol/m^2),Std Dev
    First row contains the column labels
    Subsequent rows contain data
    """
    import csv

    ## Opening the file
    textFile = open(textFileName, 'w')
    writer = csv.writer(textFile)

    writer.writerow(column_labels)
    for i in range(len(times)):
        row = [times[i], mean[i], std[i]]
        writer.writerow(row)
    textFile.close()
####

# --------------- Read ACE+ surface concentration data ----------------
def read_ACE_surf_conc(fileName):
    """Read surface concentration data stored in the standard format
    output by a reacting boundary condition in ACE+:

    # CFD-ACE+ Surface Biochemistry Output File
    # Averaged Surface Concentration (moles/m^2) for surface species
    # Patch Name : ABSORBING                                                                       
    #  Time          B(R)             AB(R)            
       0.100000E-02  0.854038501E-10  0.573649904E-12
       0.200000E-02  0.848433735E-10  0.113412651E-11
       0.300000E-02  0.843804785E-10  0.159702153E-11

    Assumes that there is only one type of adsorption site and one
    type of adsorbing species.
    """
    from numpy import empty

    f = open(fileName, 'r')
    lines = f.readlines()
    f.close()

    times = empty(len(lines)-4)
    sites = empty(len(lines)-4)
    adsorbed = empty(len(lines)-4)

    for i in range(4, len(lines)-4):
        words = lines[i].split()
        times[i] = float(words[0])
        sites[i] = float(words[1])
        adsorbed[i] = float(words[2])

    return (times, sites, adsorbed)

# --------------- Read ACE+ surface concentration data ----------------
def read_ACE_reaction_rate(fileName):
    """Read surface reaction rate data stored in the standard format
    output by a reacting boundary condition in ACE+:
    # CFD-ACE+ Surface Biochemistry Output File
    # Averaged Species Reaction Rate (moles/(m^2-sec) for all species)
    # Patch Name : RESONATORSENSINGAREA                                                            
    #  Time          A                B(R)             AB(R)            
       0.000000E+00  0.000000000E+00  0.000000000E+00  0.000000000E+00
       0.100000E-02  0.000000000E+00  0.000000000E+00  0.000000000E+00
       0.200000E-02  0.000000000E+00  0.000000000E+00  0.000000000E+00

    Assumes that there is only one type of adsorption site and one
    type of adsorbing species.
    """
    from numpy import empty

    f = open(fileName, 'r')
    lines = f.readlines()
    f.close()

    times = empty(len(lines)-4)
    bulk_species = empty(len(lines)-4)
    surface_site = empty(len(lines)-4)
    surface_complex = empty(len(lines)-4)

    for i in range(4, len(lines)-4):
        words = lines[i].split()
        times[i] = float(words[0])
        bulk_species[i] = float(words[1])
        surface_site[i] = float(words[2])
        surface_complex[i] = float(words[3])

    return (times, bulk_species, surface_site, surface_complex)

# --------------- Read ACE+ near-surface concentration data ----------------
def read_ACE_near_surface_concentration(fileName):
    """Read surface concentration data stored in the standard format
    output by a reacting boundary condition in ACE+:

    # CFD-ACE+ Surface Biochemistry Output File
    # Averaged Near Surface Concentration (mole/L) for analyte species
    # Patch Name : RESONATORSENSINGAREA                                                            
    #  Time          A
       0.000000E+00  0.000000000E+00
       0.100000E-02  0.000000000E+00
       0.200000E-02  0.000000000E+00

    Assumes that there is only one adsorbing species.
    """
    from numpy import empty, array

    f = open(fileName, 'r')
    lines = f.readlines()
    f.close()

#    times = empty(len(lines)-4)
#    conc = empty(len(lines)-4)

    times = []
    conc = []

    for i in range(4, len(lines)):
        words = lines[i].split()
#        times[i-4] = float(words[0])
#        conc[i-4] = float(words[1])
        t = float(words[0])
        c = float(words[1])
        if len(times) > 0 and t < times[-1]:
            times.pop()
            conc.pop()
        times.append(t)
        conc.append(c)

    for i in range(1, len(times)):
        if times[i-1] > times[i]:
            print 'Found t1 > t0'

    return (array(times), array(conc))

# --------------- Read ACE+ monitor file ----------------
def read_ACE_monitor_file(fileName, variable_index):
    """Read the first column of data from a monitor point file
    saved by ACE+.  Expected format (number of columns may vary):
    # TIME             A                 B(R)              AB(R)     
    #--------------------------------------------------------------------
        1.00000000E-04     0.00000000E+000   0.00000000E+000   0.00000000E+000
        2.00000000E-04     0.00000000E+000   0.00000000E+000   0.00000000E+000
        4.25462798E-04     0.00000000E+000   0.00000000E+000   0.00000000E+000
        9.36129822E-04     0.00000000E+000   0.00000000E+000   0.00000000E+000
    Arguments:
    fileName
    variable_index      selects data column (0-n)
    """
    from numpy import empty

    f = open(fileName, 'r')
    lines = f.readlines()
    f.close()

    header = lines[0].split()
    header = header[2:]     # remove # mark and TIME

    times = empty(len(lines)-2)
    conc = empty(len(lines)-2)

    for i in range(2, len(lines)):
        words = lines[i].split()
        times[i-2] = float(words[0])
        conc[i-2] = float(words[variable_index + 1])

    return (header, times, conc)

def write_2D_array(x_array, y_array, z_data, fileName):
    """Write a function of two variables z(x,y) that is sampled on a regular
    grid to a text file in a format that is easily read by a Fortran program.
    Arguments:
    x_array     1D array of x values
    y_array     1D array of y values
    z_data      2D array of data z(x,y) sampled at regular grid points
    fileName    name of file for saving data
    """
    outputFile = open(fileName, 'w')

    # x grid values
    outputFile.write(str(len(x_array)))
    outputFile.write("\n")
    for x in x_array:
        outputFile.write(str(x) + " ")
    outputFile.write("\n")

    # y grid values
    outputFile.write(str(len(y_array)))
    outputFile.write("\n")
    for y in y_array:
        outputFile.write(str(y) + " ")
    outputFile.write("\n")

    # data
    for i in range(len(x_array)):
        for j in range(len(y_array)):
            outputFile.write(str(z_data[i,j]) + " ")
        outputFile.write("\n")

    outputFile.close()

