import sys
import glob
from Bio import SeqIO

def Blat_Total_Query_Printer(name, input_psl):
    """Prints the Name and start/stops for blat matches for insertion into a genbank file"""
    List1 = Gene_Match_Combinbed_Block(input_psl)
    for entry in List1:
        if entry[0] == 0:
            entry[0] = 1
        Start = '     ' + name
        Add = 21 - len(Start)
        Line = Space_Adder(Add, Start)
        print(Line + str(entry[0]) + '..' + str(entry[1]))

def Blat_Total_Query_Lister(name, input_psl):
    """Makes a list of the Name and start/stops for blat matches for insertion into a genbank file"""
    Out = []
    List1 = Gene_Match_Combinbed_Block(input_psl)
    for entry in List1:
        if entry[0] == 0:
            entry[0] = 1
        Start = '  ' + name
        Add = 21 - len(Start)
        Line = Space_Adder(Add, Start)
        Info = Line + str(entry[0]) + '..' + str(entry[1])
        Out.append(Info)
    return Out

def Space_Adder(Length, string):
    """Adds spaces to a string to get it to a specific length"""
    Out = string
    for char in range(Length):
        Out = Out + ' '
    return Out

def Space_Before(Length, string):
    """Adds spaces before a string to get it to a specific length"""
    Out = ''
    for char in range(Length):
        Out = Out + ' '
    Out = Out + string
    return Out

def Overlap(a,b):
    return (a[0] >= b[0] and a[1] <= b[1]) or (b[0] >= a[0] and b[1] <= a[1]) or (a[0] <= b[1] and a[1] >= b[0]) or (b[0] <= a[1] and b[1] >= a[0])

##def Overlap(a,b):
##    return (a[0] >= b[0] and a[1] <= b[1]) or (b[0] >= a[0] and b[1] <= a[1]) or (a[0] == b[1] - 1) or (b[0] == a[1] - 1)


def Overlap_Extender(a, b):
    Out = a
    if (Overlap(a, b)) == True:
        if (a[0] <= b[0]):
            Out[0] = a[0]
        else:
            Out[0] = b[0]
        if a[1] <= b[1]:
            Out[1] = b[1]
        else:
            Out[1] = a[1]
    return Out

def Multiple_Overlap_Extenders(input_list):
    Out = []
##    Add = 0
    for entry in input_list:
        Focus = entry
        for entry2 in input_list:
            if Overlap(Focus, entry2) == True:
                Focus = Overlap_Extender(Focus, entry2)
##                if Focus != entry:
##                    Add = 1
        if (Focus in Out) == False:
            Out.append(Focus)
##    return [Out, Add]
    return Out

def Recursive_Overlap(input_list):
    Current = input_list
    New = Multiple_Overlap_Extenders(Current)
    while New != Current:
        Current = New
        New = Multiple_Overlap_Extenders(Current)
    return New

def Match_Start_Stop_Finder(PSL_Line):
    """Finds the start and stop for the contig and gene"""
    List1 = PSL_Line.split('\t')
    Output = []
    Block_Lengths = List1[18].split(',')[0:-1]
    Blocks = []
    for lengths in Block_Lengths:
        Blocks.append(int(lengths))
    Block_Lengths = Blocks
    Blocks = []
    Gene_Starts = List1[-1].split(',')[0:-1]
    for lengths in Gene_Starts:
        Blocks.append(int(lengths))
    Gene_Starts = Blocks
    Blocks = []
    Genome_Starts = List1[-2].split(',')[0:-1]
    for lengths in Genome_Starts:
        Blocks.append(int(lengths))
    Genome_Starts = Blocks
    if List1[8] == '-':
        Genome_Start = int(List1[10]) - int(List1[12])
        Genome_End = int(List1[10]) - int(List1[11])
    else:
        Genome_Start = int(List1[11])
        Genome_End = int(List1[12])
    Genome_Length = int(List1[10])
    Gene_Start = int(List1[15])
    Gene_End = int(List1[16])
    Gene_Length = int(List1[14])
    Gene_Starts.append(Gene_End)
    Genome_Starts.append(Genome_End)       
    Genome_Start_Stop = [Genome_Start, Genome_End]
    Output.append(Genome_Start_Stop)
    Gene_Start_Stop = [Gene_Start, Gene_End]
    Output.append(Gene_Start_Stop)
    Output.append(Block_Lengths)
    Output.append(Genome_Starts)
    Output.append(Gene_Starts)
    return Output

def Gene_Match_Block_Finder(input_psl):
    """Makes a list of all the start and stops of the overlap"""
    f = open(input_psl, 'r')
    All = []
    for line in f:
        Info = Match_Start_Stop_Finder(line)
        Blocks = Info[2]
        Starts = Info[4][0:-1]
        for entry in range(len(Blocks)):
            New = [Starts[entry], Starts[entry] + Blocks[entry]]
            All.append(New)
    return All

def Gene_Match_Combinbed_Block(input_psl):
    """Makes non-overlapping start and stop block matches"""
    Info = Gene_Match_Block_Finder(input_psl)
    Combined = Recursive_Overlap(Info)
    return Combined

def Total_Length(input_list):
    """Takes in a list of non-overlapping match blocks and returns the total length"""
    Total = 0
    for entry in input_list:
        tot = entry[1] - entry[0]
        Total += tot
    return Total

def Total_Match_Length(input_psl):
    Blocks = Gene_Match_Combinbed_Block(input_psl)
    Length = Total_Length(Blocks)
    return Length

def PSL_Plasmid_Length(input_psl):
    f = open(input_psl, 'r')
    String1 = f.readline()
    Length = int(String1.split('\t')[14])
    f.close()
    return Length

def BLAT_Overlap_Percent(input_psl):
    Match_Length = Total_Match_Length(input_psl)
    Total = PSL_Plasmid_Length(input_psl)
    Percent = Match_Length / Total
    return Percent

def Overlap_Printer(input_glob):
    List1 = glob.glob(input_glob)
    for entry in List1:
        Name = entry.split('/')[-1]
        Name = Name.split('_')[0]
        Name = Name.split('.')[0]
        Percent = BLAT_Overlap_Percent(entry)
        Length = Total_Match_Length(entry)
        Plasmid = PSL_Plasmid_Length(entry)
        print(Name + '\t' + str(Length) + '\t' + str(Percent) + '\t' + str(Plasmid))
    
