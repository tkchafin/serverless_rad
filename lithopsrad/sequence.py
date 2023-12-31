import os
import sys
import re
import io

def extract_size(header):
    """
    Extracts the size value from the fasta header.
    """
    size_tag = next((tag for tag in header.split(';') if tag.startswith('size=')), None)
    if size_tag:
        return int(size_tag.split('=')[1])
    return 0


def write_fasta(seqs, fas):
    with open(fas, 'w') as fh:
        try:
            for seq_id, sequence in seqs.items():
                line = ">" + str(seq_id) + "\n" + sequence + "\n"
                fh.write(line)
        except IOError as e:
            print("Could not read file:", e)
            sys.exit(1)
        except Exception as e:
            print("Unexpected error:", e)
            sys.exit(1)


def read_fasta(fas):
    with open(fas, "r") as file_object:
        contig = None
        seq = ""
        for line in file_object:
            line = line.strip()
            if not line:
                continue
            line = line.replace(" ","")
            #print(line)
            if line[0] == ">": #Found a header line
                #If we already loaded a contig, yield that contig and
                #start loading a new one
                if contig:
                    yield([contig,seq]) #yield
                    contig = None #reset contig and seq
                    seq = ""
                else:
                    contig = (line.replace(">",""))
                    seq = ""
            else:
                seq += line
        #yield last sequence, if it has both a header and sequence
        if contig and seq:
            yield([contig,seq])

def check_fastq(file_path):
    with open(file_path, 'r') as f:
        lines = f.readlines()
        if len(lines) % 4 != 0:
            raise ValueError(f"File {file_path} doesn't appear to be a valid FASTQ file.")


def count_fastq_from_stream(file_data):
    lines = file_data.split("\n")
    total_lines = len(lines) - lines.count('')  # Subtracting empty lines
    return total_lines // 4


def count_fastq_records(file_stream=None, file_path=None):
    if file_path is not None:
        with open(file_path, "r") as file:
            file_stream = file.read()
            return(count_fastq_from_stream(file_stream))
    else:
        return(count_fastq_from_stream(file_stream))


def get_kmers_neighbor(seq, k):
    for i in range(0, len(seq)-(k-1)):
        yield (seq[i:i+k-1], seq[i+1:i+k])


def get_kmers(seq, k):
    for i in range(0, len(seq)-(k-1)):
        yield (seq[i:i+k-1])


#Function to split character to IUPAC codes, assuing diploidy
def get_iupac_caseless(char):
    lower = False
    if char.islower():
        lower = True
        char = char.upper()
    iupac = {
        "A"    : ["A"],
        "G"    : ["G"],
        "C"    : ["C"],
        "T"    : ["T"],
        "N"    : ["A", "C", "G", "T"],
        "-"    : ["A", "C", "G", "T", "-"],
        "R"    : ["A","G"],
        "Y"    : ["C","T"],
        "S"    : ["G","C"],
        "W"    : ["A","T"],
        "K"    : ["G","T"],
        "M"    : ["A","C"],
        "B"    : ["C","G","T"],
        "D"    : ["A","G","T"],
        "H"    : ["A","C","T"],
        "V"    : ["A","C","G"]
    }
    ret = iupac[char]
    if lower:
        ret = [c.lower() for c in ret]
    return ret


#Function to expand ambiguous sequences
#Generator function
def expand_ambigs(sequence):
   for i in product(*[get_iupac_caseless(j) for j in sequence]):
       yield("".join(i))


#Function to return reverse complement of a nucleotide, while preserving case
def get_revcomp_caseless(char):
    lower = False
    if char.islower():
        lower = True
        char = char.upper()
    d = {
        "A"    : "T",
        "G"    : "C",
        "C"    : "G",
        "T"    : "A",
        "N"    : "N",
        "-"    : "-",
        "R"    : "Y",
        "Y"    : "R",
        "S"    : "S",
        "W"    : "W",
        "K"    : "M",
        "M"    : "K",
        "B"    : "V",
        "D"    : "H",
        "H"    : "D",
        "V"    : "B"
    }
    ret = d[char]
    if lower:
        ret = ret.lower()
    return ret


#Function to reverse complement a sequence, with case preserved
def revcomp(seq):
    comp = []
    for i in (get_revComp_caseless(j) for j in seq):
        comp.append(i)
    return("".join(comp[::-1]))


#Function to simplify a sequence
def simplify(seq):
    temp = re.sub('[ACGT]', '', (seq).upper())
    return temp.translate(str.maketrans("RYSWKMBDHV", "**********"))


#returns dict of character counts
def counter(seq):
    d = {}
    d = {
        'A':0,
        'N':0,
        '-':0,
        'C':0,
        'G':0,
        'T':0,
        "R"    : 0,
        "Y"    : 0,
        "S"    : 0,
        "W"    : 0,
        "K"    : 0,
        "M"    : 0,
        "B"    : 0,
        "D"    : 0,
        "H"    : 0,
        "V"    : 0
    }
    for c in seq:
        if c in d:
            d[c] += 1
    d['VAR'] = d['R'] + d['Y'] + d['S'] + d['W'] \
    + d['K'] + d['M'] + d['B'] + d['D'] + d['H'] + d['V']
    return d


#Returns dict of character counts from a simplified consensus sequence
def simple_counter(seq):
    d = {}
    d = {
        'N':0,
        '-':0,
        '*':0
    }
    for c in seq:
        if c in d:
            d[c] += 1
    return d


#Function to get GC content of a provided sequence
def gc_counts(string):
    new = re.sub('[GCgc]','#',string)
    return sum(1 for c in new if c == '#')


#Function to get counts of masked bases
def mask_counts(string):
    return sum(1 for c in string if c.islower())


#Function to get GC content as proportion
def gc_content(string):
    new = re.sub('[GCgc]','#',string)
    count = sum(1 for c in new if c == '#')
    return(count/(len(string)))


#Function to count number of lower case in a string
def mask_content(string):
    count = sum(1 for c in string if c.islower())
    return(count/(len(string)))
