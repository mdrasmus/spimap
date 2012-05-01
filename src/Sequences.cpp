/*=============================================================================

  Matt Rasmussen
  Copyright 2010-2011

  Molecular sequence input/output 

=============================================================================*/

// spidir headers
#include "seq.h"
#include "Sequences.h"
#include "parsing.h"
#include "logging.h"

namespace spidir {



Sequences *readFasta(const char *filename)
{
    FILE *infile = NULL;
    
    if ((infile = fopen(filename, "r")) == NULL) {
        printError("cannot read file '%s'", filename);
        return NULL;
    }
    
    BufferedReader reader(infile);
    char *line;
    
    Sequences *seqs = new Sequences();
    string key;
    ExtendArray<char> seq(0, 10000);

    
    while ((line = reader.readLine())) {
        chomp(line);
        
        if (line[0] == '>') {
            if (seq.size() > 0) {  
                // add new sequence
                seq.append('\0');
                seqs->append(key, seq.detach());
            }
        
            // new key found
            key = string(&line[1]);
        } else {
            seq.extend(line, strlen(line));
        }
        
    }
    
    // add last sequence
    if (seq.size() > 0) {
        seq.append('\0');
        seqs->append(key, seq.detach());
    }
    
    return seqs;
}


Sequences *readAlignFasta(const char *filename)
{
    Sequences *seq = readFasta(filename);
    if (!seq)
        return NULL;
    seq->setAlignLength();
    return seq;
}


bool writeFasta(const char *filename, Sequences *seqs)
{
    FILE *stream = NULL;
    
    if ((stream = fopen(filename, "w")) == NULL) {
        fprintf(stderr, "cannot open '%s'\n", filename);
        return false;
    }

    writeFasta(stream, seqs);
    
    fclose(stream);
    return true;
}

void writeFasta(FILE *stream, Sequences *seqs)
{
    for (int i=0; i<seqs->nseqs; i++) {
        fprintf(stream, ">%s\n", seqs->names[i].c_str());
        fprintf(stream, "%s\n", seqs->seqs[i]);
    }
}


bool checkSequences(Sequences *seqs)
{
    return checkSequences(seqs->nseqs, seqs->seqlen, seqs->seqs) &&
        checkSeqNames(seqs);
}

// ensures that all characters in the alignment are sensible
// TODO: do not change alignment (keep Ns)
bool checkSequences(int nseqs, int seqlen, char **seqs)
{
    // check seqs
    // CHANGE N's to gaps
    for (int i=0; i<nseqs; i++) {
        for (int j=0; j<seqlen; j++) {
            char x = seqs[i][j];
            if (strchr("NnRrYyWwSsKkMmBbDdHhVv", x))
                // treat Ns as gaps
                x = '-';
            if (x != '-' &&
                dna2int[(int) (unsigned char) x] == -1)
            {
                // an unknown character is in the alignment
                printError("unknown char '%c' (char code %d)\n", x, x);
                return false;
            }
        }
    }
    
    return true;
}


bool checkSeqNames(Sequences *seqs)
{
    for (int i=0; i<seqs->names.size(); i++) {
        if (!checkSeqName(seqs->names[i].c_str())) {
            printError("sequence name has illegal characters '%s'",
                       seqs->names[i].c_str());
            return false;
        }
    }

    return true;
}


//
// A valid gene name and species name follows these rules:
//
// 1. the first and last characters of the ID are a-z A-Z 0-9 _ - .
// 2. the middle characters can be a-z A-Z 0-9 _ - . or the space character ' '.
// 3. the ID should not be purely numerical characters 0-9
// 4. the ID should be unique within a gene tree or within a species tree
//

bool checkSeqName(const char *name)
{
    int len = strlen(name);
    
    if (len == 0) {
        printError("name is zero length");
        return false;
    }

    // check rule 1
    if (name[0] == ' ' || name[len-1] == ' ') {
        printError("name starts or ends with a space '%c'");
        return false;
    }

    // check rule 2
    for (int i=0; i<len; i++) {
        char c = name[i];
        if (!((c >= 'a' && c <= 'z') ||
              (c >= 'A' && c <= 'Z') ||
              (c >= '0' && c <= '9') ||
              c == '_' || c == '-' || c == '.' || c == ' ')) {
            printError("name contains illegal character '%c'", c);
            return false;
        }
    }

    // check rule 3
    int containsAlpha = false;
    for (int i=0; i<len; i++) {
        if (name[i] < '0' || name[i] > '9') {
            containsAlpha = true;
            break;
        }
    }
    if (!containsAlpha) {
        printError("name is purely numeric '%s'", name);
        return false;
    }

    return true;
}


void resampleAlign(Sequences *aln, Sequences *aln2)
{
    assert(aln->nseqs == aln2->nseqs);
    char **seqs = aln->seqs;
    char **seqs2 = aln2->seqs;

    for (int j=0; j<aln2->seqlen; j++) {
        // randomly choose a column (with replacement)
        int col = irand(aln->seqlen);
        
        // copy column
        for (int i=0; i<aln2->nseqs; i++) {
            seqs2[i][j] = seqs[i][col];
        }
    }
}




} // namespace spidir
